"""
Agent: add BIN support and richer GKA metrics for OAM imagery.

This script can ingest PNG/JPG/TIF images *and* raw .bin files, extract
ring-aligned intensity profiles, and compute GKA-like summaries
(coherence, dominant angular mode, parity, and optional BORGT slope).
"""
import argparse
from pathlib import Path
from typing import Iterable, List, Optional, Sequence, Tuple

import cv2
import numpy as np
import pandas as pd
import matplotlib
from scipy.stats import linregress
from scipy.signal import savgol_filter

# headless by default; plots are optional
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def odd_contrast_ring(intensity: np.ndarray) -> float:
    """
    Mirror-odd contrast for a 1D ring intensity I(theta).

    eta = |I(theta) - I(theta+pi)| / (0.5*(I(theta)+I(theta+pi))).
    """
    if intensity.size == 0:
        return float("nan")
    I_norm = intensity / (np.max(intensity) + 1e-9)
    shift = intensity.size // 2
    I_shift = np.roll(I_norm, shift)
    eta = np.abs(I_norm - I_shift) / (0.5 * (I_norm + I_shift) + 1e-9)
    return float(np.nanmean(eta))


def parse_shape(shape_str: Optional[str]) -> Optional[Tuple[int, ...]]:
    """Parse a shape string like '512,512' or '2,512,512' into a tuple."""
    if shape_str is None:
        return None
    cleaned = shape_str.replace("x", ",").replace(" ", "")
    parts = [p for p in cleaned.split(",") if p]
    try:
        return tuple(int(p) for p in parts)
    except ValueError:
        raise argparse.ArgumentTypeError(f"Invalid shape spec: {shape_str!r}")


def infer_square_shape(n_vals: int) -> Optional[Tuple[int, int]]:
    """Try to infer a square 2D shape from element count."""
    side = int(round(np.sqrt(n_vals)))
    if side * side == n_vals:
        return (side, side)
    return None


def load_bin_frames(
    path: Path,
    shape: Optional[Tuple[int, ...]],
    dtype: str,
    complex_mode: bool,
    max_frames: Optional[int] = None,
) -> List[np.ndarray]:
    """
    Load frames from a binary file.

    Supports:
      - real 2D arrays (ny, nx)
      - real stacks (n_frames, ny, nx)
      - complex fields (2, ny, nx) or (ny, nx, 2) when complex_mode=True
    """
    arr = np.fromfile(path, dtype=dtype)
    n_vals = arr.size

    if shape is None:
        # Try a square guess; if complex, try to peel a leading/trailing 2
        guess = infer_square_shape(n_vals)
        if guess is None and complex_mode and n_vals % 2 == 0:
            half = n_vals // 2
            sq = infer_square_shape(half)
            if sq is not None:
                shape = (2, *sq)
        else:
            shape = guess

    if shape is None:
        raise ValueError(
            f"Cannot infer shape for {path.name}; specify --bin-shape explicitly."
        )

    try:
        arr = arr.reshape(shape)
    except Exception as exc:
        raise ValueError(f"Failed to reshape {path.name} to {shape}: {exc}") from exc

    frames: List[np.ndarray] = []

    if complex_mode:
        if arr.ndim != 3:
            raise ValueError(
                f"{path.name}: complex_mode expects shape (2, ny, nx) or (ny, nx, 2); got {arr.shape}"
            )
        if arr.shape[0] == 2:
            real, imag = arr[0], arr[1]
        elif arr.shape[-1] == 2:
            real, imag = arr[..., 0], arr[..., 1]
        else:
            raise ValueError(
                f"{path.name}: cannot locate complex channels in shape {arr.shape}"
            )
        frames.append(real.astype(np.float32) ** 2 + imag.astype(np.float32) ** 2)
    elif arr.ndim == 2:
        frames.append(arr.astype(np.float32))
    elif arr.ndim == 3:
        n_available = arr.shape[0]
        limit = n_available if max_frames is None else min(max_frames, n_available)
        for idx in range(limit):
            frames.append(arr[idx].astype(np.float32))
    else:
        raise ValueError(f"{path.name}: unsupported ndim {arr.ndim}")

    return frames


def find_ring_radius(image: np.ndarray, center: Optional[Tuple[float, float]] = None) -> Tuple[float, float, float]:
    """
    Coarse ring locator using radial sampling of bright pixels.
    """
    h, w = image.shape
    if center is None:
        thresh = np.percentile(image, 90)
        ys, xs = np.where(image > thresh)
        if len(xs) == 0:
            cx, cy = w / 2.0, h / 2.0
        else:
            cx, cy = xs.mean(), ys.mean()
    else:
        cx, cy = center

    radii = np.linspace(10, min(h, w) // 2, 300)
    thetas = np.linspace(0, 2 * np.pi, 360)
    cos_t = np.cos(thetas)
    sin_t = np.sin(thetas)
    vals = []
    for r in radii:
        xs_r = (cx + r * cos_t).astype(int)
        ys_r = (cy + r * sin_t).astype(int)
        mask = (xs_r >= 0) & (xs_r < w) & (ys_r >= 0) & (ys_r < h)
        if not np.any(mask):
            vals.append(np.nan)
            continue
        vals.append(image[ys_r[mask], xs_r[mask]].mean())
    vals = np.array(vals)
    if np.all(np.isnan(vals)):
        r_best = float(min(h, w) // 4)
    else:
        r_best = float(radii[int(np.nanargmax(vals))])
    return float(cx), float(cy), float(r_best)


def extract_I_theta(image: np.ndarray, cx: float, cy: float, r: float, n_samples: int = 720) -> Tuple[np.ndarray, np.ndarray]:
    """Sample intensity along a ring at radius r with subpixel interpolation."""
    thetas = np.linspace(0, 2 * np.pi, n_samples, endpoint=False)
    xs = (cx + r * np.cos(thetas)).astype(np.float32)
    ys = (cy + r * np.sin(thetas)).astype(np.float32)
    map_x = xs.reshape(1, -1)
    map_y = ys.reshape(1, -1)
    I = cv2.remap(
        image.astype(np.float32),
        map_x,
        map_y,
        interpolation=cv2.INTER_LINEAR,
        borderMode=cv2.BORDER_REFLECT,
    ).ravel()
    return thetas, I


def quantile_contrast(vals: np.ndarray) -> float:
    """Robust contrast proxy using upper/lower quartiles."""
    p75 = np.percentile(vals, 75)
    p25 = np.percentile(vals, 25)
    return float(p75 / max(p25, 1e-6))


def gka_metrics(
    image: np.ndarray,
    cx: float,
    cy: float,
    radius: float,
    r_samples: Sequence[float],
    intensity_floor: float,
    min_radii: int,
    knee_window: int,
    save_plots: bool,
    plot_dir: Path,
    plot_label: str,
    ref_power: Optional[np.ndarray] = None,
    n_theta: int = 2048,
) -> dict:
    """Compute angular metrics for a single frame."""
    thetas, I = extract_I_theta(image, cx, cy, radius, n_samples=n_theta)
    I_norm = I / (I.max() + 1e-6)
    window = np.hanning(len(I_norm))
    I_win = I_norm * window
    median = float(np.median(I_norm))
    coherence_q = quantile_contrast(I_norm)
    coherence_max_med = float(I_norm.max() / max(median, 1e-6))

    F = np.fft.rfft(I_win)
    mags = np.abs(F)
    m_max = max(1, len(mags) // 8)
    mags_trust = mags.copy()
    mags_trust[m_max + 1 :] = 0.0
    def peak_stats(mvec: np.ndarray) -> Tuple[float, float, float]:
        if len(mvec) <= 1:
            return np.nan, np.nan, np.nan
        peak_idx = int(np.argmax(mvec[1:]) + 1)
        peak = float(mvec[peak_idx])
        second = float(np.partition(mvec[1:], -2)[-2]) if len(mvec) > 3 else 0.0
        ratio = float(peak / (second + 1e-9))
        neigh = np.delete(mvec[1:], np.argmax(mvec[1:]))
        noise = float(np.median(neigh) + 1e-9)
        z = float((peak - noise) / (np.median(np.abs(neigh - noise)) + 1e-9))
        return peak_idx, ratio, z
    m_star, m_ratio, m_z = peak_stats(mags_trust)
    energy = float((mags_trust ** 2).sum())
    purity = float((mags_trust[int(m_star)] ** 2) / (energy + 1e-9)) if np.isfinite(m_star) else np.nan
    power_odd = float((mags_trust[1::2] ** 2).sum())
    power_even = float((mags_trust[2::2] ** 2).sum())
    parity_odd = float(power_odd / (power_even + power_odd + 1e-9)) if energy > 0 else np.nan

    # Reference subtraction (intensity spectrum)
    m_star_corr = np.nan
    m_ratio_corr = np.nan
    m_z_corr = np.nan
    if ref_power is not None and len(ref_power) == len(mags_trust):
        P_img = mags_trust ** 2
        P_ref = ref_power
        fit_mask = slice(1, min(6, len(P_ref)))
        alpha = float(np.sum(P_img[fit_mask] * P_ref[fit_mask]) / (np.sum(P_ref[fit_mask] ** 2) + 1e-9))
        P_corr = P_img - alpha * P_ref
        P_corr = np.clip(P_corr, 0.0, None)
        mags_corr = np.sqrt(P_corr)
        m_star_corr, m_ratio_corr, m_z_corr = peak_stats(mags_corr)

    thr = median + 0.3 * (float(I_norm.max()) - median)
    coverage = float((I_norm > thr).mean())

    contrast_per_r = {}
    r_vals = []
    contrasts = []
    mean_norm = []
    eta_vals = []
    ring_peaks = []
    for r_here in r_samples:
        _, I_r = extract_I_theta(image, cx, cy, r_here, n_samples=n_theta)
        I_r_norm = I_r / (I_r.max() + 1e-6)
        I_r_win = I_r_norm * np.hanning(len(I_r_norm))
        c_val = quantile_contrast(I_r_norm)
        r_vals.append(r_here)
        contrasts.append(c_val)
        mean_norm.append(float(I_r_norm.mean()))
        frac = r_here / max(radius, 1e-6)
        key = f"contrast_f{int(round(frac * 100)):02d}"
        contrast_per_r[key] = c_val
        eta_vals.append(odd_contrast_ring(I_r_norm))
        mags_r = np.abs(np.fft.rfft(I_r_win))
        mags_r[m_max + 1 :] = 0.0
        pk_r, ratio_r, z_r = peak_stats(mags_r)
        ring_peaks.append((pk_r, ratio_r, z_r))

    r_vals = np.array(r_vals, dtype=np.float64)
    contrasts = np.array(contrasts, dtype=np.float64)
    mean_norm = np.array(mean_norm, dtype=np.float64)
    eta_vals = np.array(eta_vals, dtype=np.float64)

    # Ring consensus
    ring_pk_vals = [pk for pk, _, _ in ring_peaks if np.isfinite(pk)]
    ring_ratio_vals = [rt for _, rt, _ in ring_peaks if np.isfinite(rt)]
    ring_z_vals = [zt for _, _, zt in ring_peaks if np.isfinite(zt)]
    ring_consensus = float(ring_pk_vals.count(max(set(ring_pk_vals), key=ring_pk_vals.count)) / len(ring_pk_vals)) if ring_pk_vals else np.nan
    ring_ratio_med = float(np.median(ring_ratio_vals)) if ring_ratio_vals else np.nan
    ring_z_med = float(np.median(ring_z_vals)) if ring_z_vals else np.nan

    def confidence(ratio: float, z: float, consensus: float) -> float:
        if not (np.isfinite(ratio) and np.isfinite(z) and np.isfinite(consensus)):
            return 0.0
        score = 2.0 * (ratio - 1.0) + 0.5 * (z / 3.0) + 2.0 * (consensus - 0.5)
        return float(1.0 / (1.0 + np.exp(-score)))
    m_conf = confidence(m_ratio, m_z, ring_consensus)
    m_conf_corr = confidence(m_ratio_corr, m_z_corr, ring_consensus) if np.isfinite(m_ratio_corr) else 0.0

    # BORGT-like slope: log contrast vs log radius, with guards and a light outlier filter
    base_mask = (
        np.isfinite(r_vals)
        & np.isfinite(contrasts)
        & np.isfinite(mean_norm)
        & (contrasts > 0)
        & (r_vals > 0)
        & (mean_norm >= intensity_floor)
    )
    slope_status = "insufficient_radii"
    borgt_slope = np.nan
    borgt_intercept = np.nan
    slope_n = int(base_mask.sum())
    if slope_n >= min_radii:
        log_r = np.log(r_vals[base_mask])
        log_c = np.log(contrasts[base_mask])
        try:
            slope1, intercept1 = np.polyfit(log_r, log_c, 1)
            residuals = log_c - (slope1 * log_r + intercept1)
            mad = np.median(np.abs(residuals - np.median(residuals))) + 1e-9
            keep = np.abs(residuals) <= 2.5 * mad
            slope_n = int(keep.sum())
            if slope_n >= min_radii:
                slope2, intercept2 = np.polyfit(log_r[keep], log_c[keep], 1)
                borgt_slope = float(slope2)
                borgt_intercept = float(intercept2)
                slope_status = "ok"
            else:
                slope_status = "outlier_filtered_all"
        except Exception:
            slope_status = "fit_failed"
    elif slope_n > 0:
        slope_status = "low_signal"

    # Mirror-odd slope vs radius
    eta_mask = (
        np.isfinite(r_vals)
        & np.isfinite(eta_vals)
        & np.isfinite(mean_norm)
        & (eta_vals > 0)
        & (r_vals > 0)
        & (mean_norm >= intensity_floor)
    )
    eta_status = "insufficient_radii"
    eta_slope = np.nan
    eta_intercept = np.nan
    eta_r2 = np.nan
    eta_knee = np.nan
    eta_n = int(eta_mask.sum())
    if eta_n >= min_radii:
        log_r_eta = np.log(r_vals[eta_mask])
        log_eta = np.log(eta_vals[eta_mask])
        try:
            sl, itc, r_val, _, se_val = linregress(log_r_eta, log_eta)
        except Exception:
            eta_status = "fit_failed"
        else:
            eta_slope = float(sl)
            eta_intercept = float(itc)
            eta_r2 = float(r_val**2)
            eta_status = "ok"
            # knee via smoothed curvature on log-log
            size_ok = log_r_eta.size
            if size_ok >= max(5, knee_window):
                try:
                    w = knee_window if knee_window % 2 == 1 else knee_window + 1
                    if w >= size_ok:
                        w = size_ok if size_ok % 2 == 1 else size_ok - 1
                    d1 = np.gradient(log_eta, log_r_eta)
                    d1s = savgol_filter(d1, window_length=max(3, w), polyorder=2)
                    d2 = np.gradient(d1s, log_r_eta)
                    idx = int(np.argmax(np.abs(d2)))
                    eta_knee = float(r_vals[eta_mask][idx])
                except Exception:
                    eta_knee = np.nan
            else:
                eta_status = "knee_insufficient"
    elif eta_n > 0:
        eta_status = "low_signal"

    if save_plots:
        plot_dir.mkdir(parents=True, exist_ok=True)
        plt.figure(figsize=(8, 3))
        plt.plot(thetas, I_norm)
        plt.title(f"I(theta) | m*={m_star}")
        plt.xlabel("theta (rad)")
        plt.ylabel("Intensity (norm)")
        plt.tight_layout()
        plt.savefig(plot_dir / f"{plot_label}_itheta.png", dpi=150)
        plt.close()

    return {
        "cx": cx,
        "cy": cy,
        "radius": radius,
        "coherence_q": coherence_q,
        "coherence_max_med": coherence_max_med,
        "m_star": m_star,
        "m_ratio": m_ratio,
        "m_z": m_z,
        "m_star_corr": m_star_corr,
        "m_ratio_corr": m_ratio_corr,
        "m_z_corr": m_z_corr,
        "m_conf": m_conf,
        "m_conf_corr": m_conf_corr,
        "ring_consensus": ring_consensus,
        "ring_ratio_med": ring_ratio_med,
        "ring_z_med": ring_z_med,
        "purity": purity,
        "parity_odd": parity_odd,
        "coverage": coverage,
        "borgt_slope": borgt_slope,
        "borgt_intercept": borgt_intercept,
        "borgt_status": slope_status,
        "borgt_n": slope_n,
        "eta_slope": eta_slope,
        "eta_intercept": eta_intercept,
        "eta_r2": eta_r2,
        "eta_status": eta_status,
        "eta_n": eta_n,
        "eta_knee": eta_knee,
        **contrast_per_r,
    }


def normalize_image(img: np.ndarray) -> np.ndarray:
    """Normalize to float32 [0, 1]."""
    img_f = img.astype(np.float32)
    max_val = img_f.max()
    if max_val > 0:
        img_f = img_f / max_val
    return img_f


def collect_images(root: Path, pattern: str) -> List[Path]:
    """Recursively collect image paths matching the pattern."""
    return sorted(root.rglob(pattern))


def main():
    parser = argparse.ArgumentParser(
        description="GKA/OAM analysis with image + BIN support (quantile coherence, BORGT slope, parity)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--root", default="Data and code", help="root directory to scan")
    parser.add_argument("--image-glob", default="*.png", help="glob for images (relative to root, recursive)")
    parser.add_argument("--bin-pattern", default=None, help="glob for .bin files (relative to root, recursive)")
    parser.add_argument("--bin-shape", type=parse_shape, default=None, help="explicit shape for .bin (e.g. '1024,1024' or '2,1024,1024')")
    parser.add_argument("--bin-dtype", default="uint16", help="dtype for .bin files (e.g. uint16, float32)")
    parser.add_argument("--bin-complex", action="store_true", help="interpret .bin as complex field (Re/Im)")
    parser.add_argument("--bin-max-frames", type=int, default=None, help="limit frames per bin file (if stack)")
    parser.add_argument("--r-min", type=float, default=0.5, help="minimum radius fraction (relative to detected ring) for radial sampling")
    parser.add_argument("--r-max", type=float, default=0.95, help="maximum radius fraction for radial sampling")
    parser.add_argument("--r-count", type=int, default=9, help="number of radii to sample between r-min and r-max")
    parser.add_argument("--intensity-floor", type=float, default=0.02, help="minimum mean ring intensity (normalized) to keep a radius in slope fit")
    parser.add_argument("--min-radii", type=int, default=6, help="minimum # of good radii required for BORGT slope fit")
    parser.add_argument("--knee-window", type=int, default=7, help="window (odd) for knee smoothing on log-log curves")
    parser.add_argument("--out-csv", default="gka_oam_image_summary_v4.csv", help="output CSV path")
    parser.add_argument("--plots", action="store_true", help="save I(theta) plots")
    parser.add_argument("--plot-dir", default=None, help="optional plot directory (defaults to <root>/gka_oam_plots)")
    parser.add_argument("--n-theta", type=int, default=2048, help="samples along theta for ring extraction")

    args = parser.parse_args()
    root = Path(args.root)

    if args.r_count <= 1:
        r_fractions = [args.r_min]
    else:
        r_fractions = np.linspace(args.r_min, args.r_max, args.r_count)

    plot_dir = Path(args.plot_dir) if args.plot_dir else (root / "gka_oam_plots")

    records = []

    # Images
    img_paths = collect_images(root, args.image_glob)
    ref_power_cache: dict = {}
    # First pass: collect reference spectra (filenames containing inphase/zero/ref)
    keywords = ("inphase", "reference", "ref_", "zero")
    for path in img_paths:
        lower = path.name.lower()
        if any(k in lower for k in keywords):
            img = cv2.imread(str(path), cv2.IMREAD_GRAYSCALE)
            if img is None:
                continue
            img_norm = normalize_image(img)
            cx, cy, r = find_ring_radius(img_norm)
            _, I_ref = extract_I_theta(img_norm, cx, cy, r, n_samples=args.n_theta)
            I_ref_norm = I_ref / (I_ref.max() + 1e-6)
            I_ref_win = I_ref_norm * np.hanning(len(I_ref_norm))
            m_ref = np.abs(np.fft.rfft(I_ref_win))
            m_max_ref = max(1, len(m_ref) // 8)
            m_ref[m_max_ref + 1 :] = 0.0
            ref_power_cache[path.parent] = (m_ref ** 2)

    for path in img_paths:
        img = cv2.imread(str(path), cv2.IMREAD_GRAYSCALE)
        if img is None:
            continue
        img_norm = normalize_image(img)
        cx, cy, r = find_ring_radius(img_norm)
        ref_power = ref_power_cache.get(path.parent)
        metrics = gka_metrics(
            img_norm,
            cx,
            cy,
            r,
            r * r_fractions,
            args.intensity_floor,
            args.min_radii,
            args.knee_window,
            args.plots,
            plot_dir,
            plot_label=path.stem,
            ref_power=ref_power,
            n_theta=args.n_theta,
        )
        records.append(
            {
                "source": str(path.relative_to(root)),
                "kind": "image",
                **metrics,
            }
        )

    # BIN files
    if args.bin_pattern:
        bin_paths = collect_images(root, args.bin_pattern)
        for path in bin_paths:
            try:
                frames = load_bin_frames(
                    path,
                    shape=args.bin_shape,
                    dtype=args.bin_dtype,
                    complex_mode=args.bin_complex,
                    max_frames=args.bin_max_frames,
                )
            except Exception as exc:
                print(f"Skipping {path} (bin load error): {exc}")
                continue

            for idx, frame in enumerate(frames):
                frame_norm = normalize_image(frame)
                cx, cy, r = find_ring_radius(frame_norm)
                metrics = gka_metrics(
                    frame_norm,
                    cx,
                    cy,
                    r,
                    r * r_fractions,
                    args.intensity_floor,
                    args.min_radii,
                    args.knee_window,
                    args.plots,
                    plot_dir,
                    plot_label=f"{path.stem}_{idx}",
                    ref_power=None,
                    n_theta=args.n_theta,
                )
                records.append(
                    {
                        "source": f"{path.relative_to(root)}[{idx}]",
                        "kind": "bin",
                        **metrics,
                    }
                )

    if not records:
        raise SystemExit("No inputs processed; check your glob patterns.")

    df = pd.DataFrame(records)
    out_path = Path(args.out_csv)
    if not out_path.is_absolute():
        out_path = root / out_path
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)

    print(f"Processed {len(records)} items. Summary -> {out_path}")
    print(df.head())


if __name__ == "__main__":
    main()
