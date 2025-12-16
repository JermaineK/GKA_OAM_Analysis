"""
Lightweight, input-agnostic spiral/OAM witness runner for real images.

This script loads arbitrary images/arrays, estimates a center, builds a
phase-less complex surrogate field, and reports helical spectrum stability
plus intensity-spectrum diagnostics. Designed for Talbot/OAM PNGs but works
for any grayscale input.
"""
import re
from pathlib import Path
from typing import Dict, List, Tuple

import argparse
import hashlib
import random
from pathlib import Path
from typing import List, Tuple

import cv2
import numpy as np

from spiral_decoder_experiment import (
    estimate_center,
    helical_spectrum_stability,
    detect_sideband,
    reconstruct_complex_field,
)


def load_any(path: Path) -> np.ndarray:
    """Try to load grayscale image; fall back to numpy array."""
    img = cv2.imread(str(path), cv2.IMREAD_GRAYSCALE)
    if img is not None:
        return img.astype(np.float32)
    # numpy fallback
    arr = np.load(path)
    if isinstance(arr, np.lib.npyio.NpzFile):
        # pick first array
        for k in arr.files:
            return arr[k].astype(np.float32)
        raise ValueError(f"{path}: npz is empty")
    return arr.astype(np.float32)


def normalize(img: np.ndarray) -> np.ndarray:
    mx = float(np.max(img))
    if mx <= 0:
        return np.zeros_like(img, dtype=np.float32)
    return (img / mx).astype(np.float32)


def parse_condition_labels(path: Path) -> Dict[str, float]:
    """Parse weak labels like '12-to-2' from any part of the path."""
    txt = str(path).replace("\\", "/")
    m = re.search(r"(-?\d+)[-_]*to[_-]*(-?\d+)", txt, flags=re.IGNORECASE)
    if not m:
        return {
            "label_l1": np.nan,
            "label_l2": np.nan,
            "label_delta": np.nan,
            "label_sum_abs": np.nan,
            "label_abs_l1": np.nan,
            "label_abs_l2": np.nan,
        }
    l1 = float(m.group(1))
    l2 = float(m.group(2))
    return {
        "label_l1": l1,
        "label_l2": l2,
        "label_delta": abs(l1 - l2),
        "label_sum_abs": abs(l1 + l2),
        "label_abs_l1": abs(l1),
        "label_abs_l2": abs(l2),
    }


def intensity_spectrum(
    img: np.ndarray, center: Tuple[float, float], r_frac: float, n_theta: int = 2048, return_signal: bool = False
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute intensity I(theta) spectrum with angular window and oversampling at a given radius fraction.

    If return_signal is True, also return the normalized ring samples and window used.
    """
    h, w = img.shape
    cx, cy = center
    r = r_frac * (min(h, w) / 2.0)
    thetas = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
    xs = (cx + r * np.cos(thetas)).astype(np.float32)
    ys = (cy + r * np.sin(thetas)).astype(np.float32)
    I = cv2.remap(img, xs.reshape(1, -1), ys.reshape(1, -1), interpolation=cv2.INTER_LINEAR, borderMode=cv2.BORDER_REFLECT).ravel()
    I = I / (I.max() + 1e-6)
    win = np.hanning(len(I))
    Iw = I * win
    mags = np.abs(np.fft.rfft(Iw))
    m_max = max(1, len(mags) // 8)
    mags[m_max + 1 :] = 0.0
    m_axis = np.arange(mags.size)
    if return_signal:
        return m_axis, mags, I, win
    return m_axis, mags


def peak_stats(mags: np.ndarray) -> Tuple[float, float, float]:
    """Return peak m, peak ratio, and spectral entropy."""
    if mags.size <= 1:
        return np.nan, np.nan, np.nan
    peak_idx = int(np.argmax(mags[1:]) + 1)
    vals = mags[1:].astype(np.float64)
    if vals.size > 1:
        second = float(np.partition(vals, -2)[-2])
    else:
        second = 0.0
    ratio = float(vals[peak_idx - 1] / (second + 1e-9))
    p = vals / (np.sum(vals) + 1e-12)
    entropy = float(-(p * np.log(p + 1e-12)).sum())
    return float(peak_idx), ratio, entropy


def abs_spectrum_from_spec(spec: dict) -> Tuple[float, float, float]:
    """Collapse signed ell entries into |ell| groups; return peak |ell|, ratio, entropy."""
    ell_keys = [k for k in spec.keys() if k.startswith("spec_med_ell_")]
    abs_map = {}
    for k in ell_keys:
        ell = int(k.split("spec_med_ell_")[-1])
        val = float(spec.get(k, 0.0))
        abs_map.setdefault(abs(ell), 0.0)
        abs_map[abs(ell)] += val
    if not abs_map:
        return np.nan, np.nan, np.nan
    items = sorted(abs_map.items(), key=lambda t: t[1], reverse=True)
    peak_abs, peak_val = items[0]
    second_val = items[1][1] if len(items) > 1 else 1e-9
    ratio = float(peak_val / (second_val + 1e-9))
    vals = np.array([v for _, v in items], dtype=np.float64)
    p = vals / (np.sum(vals) + 1e-12)
    entropy = float(-(p * np.log(p + 1e-12)).sum())
    return float(peak_abs), ratio, entropy


def sign_observable(spec: dict, tol: float = 1e-3) -> bool:
    """Check if spec shows meaningful asymmetry between +ell and -ell."""
    ell_keys = [k for k in spec.keys() if k.startswith("spec_med_ell_")]
    if not ell_keys:
        return False
    ells = sorted({int(k.split("spec_med_ell_")[-1]) for k in ell_keys})
    for ell in ells:
        if ell <= 0:
            continue
        pos = float(spec.get(f"spec_med_ell_{ell}", 0.0))
        neg = float(spec.get(f"spec_med_ell_{-ell}", 0.0))
        if abs(pos - neg) > tol:
            return True
    return False


def helical_coeffs(field: np.ndarray, center: Tuple[float, float], ells: Tuple[int, ...], r_inner_frac: float, r_outer_frac: float, offsets: Tuple[Tuple[float, float], ...]) -> dict:
    """Compute complex helical projection coefficients per ell across jitters."""
    h, w = field.shape
    yy, xx = np.mgrid[0:h, 0:w]
    results = {ell: [] for ell in ells}
    for dx, dy in offsets:
        cx = center[0] + dx
        cy = center[1] + dy
        r = np.hypot(xx - cx, cy - yy)
        r_min = r_inner_frac * (min(h, w) / 2.0)
        r_max = r_outer_frac * (min(h, w) / 2.0)
        ann = (r >= r_min) & (r <= r_max)
        if not np.any(ann):
            continue
        th = np.arctan2(-(yy - cy), (xx - cx))
        z = field[ann]
        a = np.abs(z)
        norm = float(np.sum(a) + 1e-9)
        for ell in ells:
            basis = np.exp(-1j * ell * th[ann])
            c = np.sum(z * basis) / norm
            results[ell].append(c)
    return results


def phase_slope_signed(field: np.ndarray, center: Tuple[float, float], r_inner: float, r_outer: float, offsets: Tuple[Tuple[float, float], ...], n_theta: int = 512) -> dict:
    """Estimate signed ell via phase slope around a ring; median across offsets."""
    h, w = field.shape
    slopes = []
    r_mid = 0.5 * (r_inner + r_outer) * (min(h, w) / 2.0)
    thetas = np.linspace(-np.pi, np.pi, n_theta, endpoint=False)
    win = np.hanning(n_theta)
    for dx, dy in offsets:
        cx = center[0] + dx
        cy = center[1] + dy
        xs = (cx + r_mid * np.cos(thetas)).astype(np.float32)
        ys = (cy - r_mid * np.sin(thetas)).astype(np.float32)
        real = cv2.remap(np.real(field).astype(np.float32), xs.reshape(1, -1), ys.reshape(1, -1), cv2.INTER_LINEAR, borderMode=cv2.BORDER_REFLECT).ravel()
        imag = cv2.remap(np.imag(field).astype(np.float32), xs.reshape(1, -1), ys.reshape(1, -1), cv2.INTER_LINEAR, borderMode=cv2.BORDER_REFLECT).ravel()
        u = real + 1j * imag
        amp = np.abs(u)
        if np.median(amp) <= 1e-6:
            continue
        phi = np.unwrap(np.angle(u))
        # weighted linear fit phi ~ k * theta
        wgt = amp * win
        wgt = wgt / (np.sum(wgt) + 1e-9)
        t0 = np.sum(thetas * wgt)
        p0 = np.sum(phi * wgt)
        num = np.sum(wgt * (thetas - t0) * (phi - p0))
        den = np.sum(wgt * (thetas - t0) ** 2) + 1e-12
        slope = num / den
        ell_hat = slope  # radians per rad
        slopes.append(ell_hat)
    if not slopes:
        return {"phase_ell_med": float("nan"), "phase_ell_std": float("nan"), "phase_ell_sign_consensus": float("nan")}
    slopes = np.array(slopes, dtype=np.float64)
    signs = np.sign(slopes)
    sign_consensus = float(np.sum(signs == np.sign(np.median(slopes))) / len(signs))
    return {
        "phase_ell_med": float(np.median(slopes)),
        "phase_ell_std": float(np.std(slopes)),
        "phase_ell_sign_consensus": sign_consensus,
    }


def reconstruct_offaxis_sweep(img: np.ndarray) -> Tuple[np.ndarray, Dict[str, float]]:
    """Run a small sweep over carrier_radius_frac to reduce twin leak; returns best field and diagnostics."""
    h, w = img.shape
    # background suppression to reduce DC
    img_u8 = np.clip(img * 255.0 / (img.max() + 1e-6), 0, 255).astype(np.uint8)
    bg_u8 = cv2.medianBlur(img_u8, ksize=31)
    bg = bg_u8.astype(np.float32) / 255.0
    pre = img - bg
    pre = pre - pre.min()
    pre = pre / (pre.max() + 1e-6)
    window = np.outer(np.hanning(h), np.hanning(w)).astype(np.float32)
    F = np.fft.fftshift(np.fft.fft2(pre * window))
    (py_det, px_det), radius_det, peak_snr_det, _ = detect_sideband(np.abs(F))
    cy_fft, cx_fft = h // 2, w // 2
    fx_pred = (px_det - cx_fft) / w
    fy_pred = (py_det - cy_fft) / h
    ref_kx = float(fx_pred * 2 * np.pi)
    ref_ky = float(fy_pred * 2 * np.pi)
    best = None
    best_diag = {}
    candidates = []
    for r_frac in (0.08, 0.12, 0.18, 0.24, 0.32):
        try:
            comp, diag = reconstruct_complex_field(pre, ref_kx=ref_kx, ref_ky=ref_ky, carrier_radius_frac=r_frac, auto_sideband=True)
            diag = {**diag, "used_r_frac": r_frac}
            candidates.append((diag.get("twin_leak", 1.0), -diag.get("sideband_snr", 0.0), comp, diag))
        except Exception:
            continue
    if candidates:
        candidates.sort(key=lambda t: (t[0], t[1]))
        best = candidates[0][2]
        best_diag = candidates[0][3]
    if best is None:
        # fallback: single reconstruct
        comp, diag = reconstruct_complex_field(pre, ref_kx=ref_kx, ref_ky=ref_ky, auto_sideband=True)
        best = comp
        best_diag = {**diag, "used_r_frac": np.nan}
    best_diag.update(
        {
            "det_peak_x": float(px_det),
            "det_peak_y": float(py_det),
            "det_radius": float(radius_det),
            "det_snr": float(peak_snr_det),
            "ref_kx_hat": ref_kx,
            "ref_ky_hat": ref_ky,
        }
    )
    return best.astype(np.complex64), best_diag


def process_path(path: Path, out_dir: Path, n_theta: int, rings: Tuple[float, float], field_kind: str) -> dict:
    raw = load_any(path)
    img = normalize(raw)
    center = estimate_center(img)
    labels = parse_condition_labels(path)
    ells = tuple(range(-32, 0)) + tuple(range(1, 33))
    offsets = ((0, 0), (1, 0), (-1, 0), (0, 1), (0, -1))

    # Attempt complex reconstruction when requested
    field = img.astype(np.complex64)
    recon_diag = {}
    field_complex = field_kind in {"complex", "demod"}
    if field_kind == "complex":
        field = raw.astype(np.complex64)
    elif field_kind == "demod":
        try:
            comp, diag = reconstruct_offaxis_sweep(img)
            field = comp.astype(np.complex64)
            recon_diag = {f"recon_{k}": v for k, v in diag.items() if isinstance(v, (float, int))}
        except Exception as exc:
            recon_diag = {"recon_error": str(exc)}
            field_complex = False

    # helical spectrum stability on phase-only unit phasor
    spec = helical_spectrum_stability(
        field=field / (np.abs(field) + 1e-9),
        center=center,
        ells=ells,
        bands=((rings[0], rings[1]),),
        center_offsets=offsets,
    )
    # whiten/z-score spectrum per image to reduce persistent biases
    med_vals = np.array([spec[f"spec_med_ell_{ell}"] for ell in ells], dtype=np.float64)
    med_all = float(np.median(med_vals))
    mad_all = float(np.median(np.abs(med_vals - med_all)) + 1e-9)
    resid_vals = {ell: (spec[f"spec_med_ell_{ell}"] - med_all) / mad_all for ell in ells}
    resid_items = sorted(resid_vals.items(), key=lambda t: t[1], reverse=True)
    resid_peak_ell, resid_peak_val = resid_items[0]
    resid_second_val = resid_items[1][1] if len(resid_items) > 1 else 1e-9
    resid_peak_ratio = float(resid_peak_val / (resid_second_val + 1e-9))
    spec_mean = float(np.mean(med_vals))
    spec_std = float(np.std(med_vals) + 1e-9)
    spec_entropy = float(-(np.array(med_vals) / (np.sum(med_vals) + 1e-12) * np.log((np.array(med_vals) / (np.sum(med_vals) + 1e-12)) + 1e-12)).sum())
    spec_peak_z = float((spec.get("spec_peak_val", 0.0) - spec_mean) / spec_std)
    abs_peak, abs_ratio, abs_entropy = abs_spectrum_from_spec(spec)
    coeffs = {}
    top_modes = []
    sign_pref_peak = float("nan")
    sign_stability_peak = float("nan")
    pair_unsigned = {}
    pair_signed = {}
    pair_signed_stab = {}
    sign_ok = False
    phase_diag = {"phase_ell_med": float("nan"), "phase_ell_std": float("nan"), "phase_ell_sign_consensus": float("nan")}
    if field_complex:
        coeffs = helical_coeffs(
            field=field,
            center=center,
            ells=ells,
            r_inner_frac=rings[0],
            r_outer_frac=rings[1],
            offsets=offsets,
        )
        # top-K |c_ell| ordering
        mags = []
        for ell, vals in coeffs.items():
            if not vals:
                continue
            mags.append((ell, float(np.median(np.abs(vals)))))
        mags_sorted = sorted(mags, key=lambda t: t[1], reverse=True)
        top_modes = mags_sorted[:3]
        # pair unsigned/signed summaries
        for k in range(1, 13):
            pos = np.array([c for c in coeffs.get(k, [])], dtype=np.complex128)
            neg = np.array([c for c in coeffs.get(-k, [])], dtype=np.complex128)
            if pos.size == 0 and neg.size == 0:
                continue
            u = float(np.median(np.abs(pos))) + float(np.median(np.abs(neg)))
            pair_unsigned[k] = u
            if pos.size and neg.size:
                mlen = min(len(pos), len(neg))
                prefs = np.abs(pos[:mlen]) - np.abs(neg[:mlen])
                pair_signed[k] = float(np.median(prefs))
                pair_signed_stab[k] = float(np.mean(prefs > 0))
        if pair_unsigned:
            peak_abs_k = max(pair_unsigned.items(), key=lambda t: t[1])[0]
            sign_pref_peak = pair_signed.get(peak_abs_k, float("nan"))
            sign_stability_peak = pair_signed_stab.get(peak_abs_k, float("nan"))
            asym_max = max(abs(v) for v in pair_signed.values()) if pair_signed else 0.0
            sign_ok = bool(asym_max > 1e-3)
        phase_diag = phase_slope_signed(
            field=field,
            center=center,
            r_inner=rings[0],
            r_outer=rings[1],
            offsets=offsets,
            n_theta=512,
        )

    # intensity spectrum at nominal ring (harmonic channel)
    r_nom = rings[1]
    m_axis, mags, ring_signal, win = intensity_spectrum(img, center=center, r_frac=r_nom, n_theta=n_theta, return_signal=True)
    m_peak, m_ratio, m_entropy = peak_stats(mags)
    # bandpower split
    low_band = float(np.sum(mags[:7] ** 2))  # |m|<=6
    high_band = float(np.sum(mags[7:] ** 2))
    # ring width sweep stability
    r_samples = np.linspace(rings[0], rings[1], 4)
    ring_peaks = []
    ring_ratios = []
    for r_frac in r_samples:
        _, mags_r = intensity_spectrum(img, center=center, r_frac=float(r_frac), n_theta=n_theta)
        pk_r, ratio_r, _ = peak_stats(mags_r)
        ring_peaks.append(pk_r)
        ring_ratios.append(ratio_r)
    ring_consensus = float(ring_peaks.count(max(set(ring_peaks), key=ring_peaks.count)) / len(ring_peaks)) if ring_peaks else np.nan
    ring_ratio_med = float(np.median(ring_ratios)) if ring_ratios else np.nan
    # center jitter sensitivity
    jitter_peaks = []
    for dx, dy in offsets:
        c_j = (center[0] + dx, center[1] + dy)
        _, mags_j = intensity_spectrum(img, center=c_j, r_frac=r_nom, n_theta=n_theta)
        pk_j, _, _ = peak_stats(mags_j)
        jitter_peaks.append(pk_j)
    jitter_consensus = float(jitter_peaks.count(max(set(jitter_peaks), key=jitter_peaks.count)) / len(jitter_peaks)) if jitter_peaks else np.nan
    # angular shuffle null (destroys azimuthal coherence)
    obs_peak_val = float(np.max(mags[1:])) if mags.size > 1 else float("nan")
    null_peaks = []
    for _ in range(20):
        shuffled = np.random.permutation(ring_signal)
        mags_n = np.abs(np.fft.rfft(shuffled * win))
        m_max = max(1, len(mags_n) // 8)
        mags_n[m_max + 1 :] = 0.0
        if mags_n.size > 1:
            null_peaks.append(float(np.max(mags_n[1:])))
    null_p = float(np.mean(np.array(null_peaks) >= obs_peak_val)) if null_peaks else float("nan")
    # alias warnings
    alias_warn_spec = bool(abs(float(spec.get("spec_peak_ell", 0.0))) >= 0.9 * max(abs(ells)))
    alias_warn_harm = bool(m_peak >= (len(m_axis) - 1) * 0.9)
    # confidence tiers
    def tier(ratio: float, z: float) -> str:
        if ratio >= 1.35 and z >= 3.0:
            return "confident"
        if ratio >= 1.20 and z >= 2.0:
            return "probable"
        return "ambiguous"
    spec_confidence = tier(spec.get("spec_peak_ratio", 0.0), spec_peak_z)
    harm_mean = float(np.mean(mags[1:])) if mags.size > 1 else 0.0
    harm_std = float(np.std(mags[1:]) + 1e-9) if mags.size > 1 else 1e-9
    harm_peak_val = float(mags[int(m_peak)]) if m_peak == m_peak else 0.0
    harm_peak_z = float((harm_peak_val - harm_mean) / harm_std)
    harm_confidence = tier(m_ratio, harm_peak_z)
    # save diagnostic plots/images
    out_dir.mkdir(parents=True, exist_ok=True)
    hsh = hashlib.md5(str(path).encode("utf-8")).hexdigest()[:8]
    norm_png = out_dir / f"{path.stem}_{hsh}_norm.png"
    cv2.imwrite(str(norm_png), (img * 255).astype(np.uint8))
    spec_png = out_dir / f"{path.stem}_{hsh}_spectrum.png"
    try:
        import matplotlib.pyplot as plt

        plt.figure(figsize=(6, 3))
        plt.plot(m_axis, mags)
        plt.title(f"I(theta) spectrum | {path.name}")
        plt.xlabel("m")
        plt.ylabel("|F|")
        plt.tight_layout()
        plt.savefig(spec_png, dpi=150)
        plt.close()
    except Exception:
        pass
    recon_twin = recon_diag.get("recon_twin_leak", np.inf)
    quality_pass = bool(
        (recon_twin < 0.6)
        and (spec.get("spec_stability_frac", 0.0) > 0.85)
        and (spec.get("spec_peak_ratio", 0.0) > 1.25)
    )
    return {
        "source": str(path),
        "cx": center[0],
        "cy": center[1],
        **spec,
        "spec_resid_peak_ell": float(resid_peak_ell),
        "spec_resid_peak_val": float(resid_peak_val),
        "spec_resid_peak_ratio": float(resid_peak_ratio),
        "abs_peak_ell": abs_peak,
        "abs_peak_ratio": abs_ratio,
        "abs_entropy": abs_entropy,
        "sign_observable": bool(sign_ok),
        "sign_pref_peak": sign_pref_peak,
        "sign_stability_peak": sign_stability_peak,
        "pair_unsigned": str(pair_unsigned),
        "pair_signed": str(pair_signed),
        "pair_signed_stab": str(pair_signed_stab),
        "int_peak_m": m_peak,
        "int_peak_ratio": m_ratio,
        "int_entropy": m_entropy,
        "bandpower_low": low_band,
        "bandpower_high": high_band,
        "ring_peak_consensus": ring_consensus,
        "ring_ratio_median": ring_ratio_med,
        "center_jitter_consensus": jitter_consensus,
        "null_p_shuffle": null_p,
        "top_modes": str(top_modes),
        **phase_diag,
        "field_kind": field_kind,
        **recon_diag,
        "quality_pass": quality_pass,
        **labels,
        "spec_plot": str(spec_png),
        "norm_image": str(norm_png),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description="Real-image spiral/OAM witness runner (agnostic inputs).")
    ap.add_argument("--inputs", nargs="+", help="Explicit file paths.")
    ap.add_argument("--root", type=Path, default=Path("Generalized angle-OAM Talbot effect"), help="Root directory to scan when --glob is used.")
    ap.add_argument("--glob", default="**/*.png", help="Glob pattern under root.")
    ap.add_argument("--out-dir", type=Path, default=Path("Results/real_witness"), help="Output directory for diagnostics.")
    ap.add_argument("--n-theta", type=int, default=2048, help="Angular samples for spectra.")
    ap.add_argument("--r-inner", type=float, default=0.05, help="Inner radius fraction for annulus.")
    ap.add_argument("--r-outer", type=float, default=0.45, help="Outer radius fraction for annulus.")
    ap.add_argument(
        "--field-kind",
        choices=["intensity", "complex", "demod"],
        default="intensity",
        help="Interpret input as intensity-only, provided complex field, or demodulated from off-axis carrier.",
    )
    args = ap.parse_args()

    paths: List[Path] = []
    if args.inputs:
        paths.extend([Path(p) for p in args.inputs])
    else:
        paths.extend(sorted(args.root.glob(args.glob)))
    paths = [p for p in paths if p.exists()]
    if not paths:
        raise SystemExit("No inputs found.")

    rows = []
    for p in paths:
        try:
            row = process_path(p, out_dir=args.out_dir, n_theta=args.n_theta, rings=(args.r_inner, args.r_outer), field_kind=args.field_kind)
            rows.append(row)
        except Exception as exc:
            print(f"Failed on {p}: {exc}")

    if not rows:
        raise SystemExit("No results produced.")
    import pandas as pd

    df = pd.DataFrame(rows)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.out_dir / "real_witness_results.csv"
    df.to_csv(csv_path, index=False)
    print(f"Wrote {len(df)} rows to {csv_path}")


if __name__ == "__main__":
    main()
