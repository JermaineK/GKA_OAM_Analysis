"""
Apply symmetry-driven filters to image datasets and produce real vs null baselines.

For each image in given roots:
- detect ring radius (reuse gka_oam_extended.find_ring_radius)
- sample multiple radii (default 5) to build I(r, theta)
- compute parity-odd contrast and spiral-harmonic contrast (symmetry_filters)
- compute knee on log-log for eta_odd
- compute same metrics on a null (theta-shuffled) version

Outputs: Results/summaries/symmetry_meta.csv
"""
from pathlib import Path
from typing import Sequence
import numpy as np
import pandas as pd
import cv2

from gka_oam_extended import find_ring_radius, extract_I_theta, normalize_image
import symmetry_filters as sf


def sample_I_r_theta(img: np.ndarray, cx: float, cy: float, r_base: float, r_fractions: Sequence[float], n_theta: int = 360) -> (np.ndarray, np.ndarray):
    """Build I(r,theta) for multiple radii."""
    radii = r_base * np.asarray(r_fractions)
    thetas = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
    I_rows = []
    for r in radii:
        _, I = extract_I_theta(img, cx, cy, r)
        I_rows.append(I)
    return np.vstack(I_rows), radii


def theta_shuffle(I_r_theta: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Shuffle theta dimension (independently per radius)."""
    out = np.empty_like(I_r_theta)
    for i in range(I_r_theta.shape[0]):
        perm = rng.permutation(I_r_theta.shape[1])
        out[i] = I_r_theta[i, perm]
    return out


def process_image(path: Path, r_fractions, rng) -> dict:
    img = cv2.imread(str(path), cv2.IMREAD_GRAYSCALE)
    if img is None:
        raise ValueError(f"Failed to read {path}")
    # downsample large images for speed
    h, w = img.shape
    max_dim = max(h, w)
    if max_dim > 1024:
        scale = 1024.0 / max_dim
        img = cv2.resize(img, (int(w * scale), int(h * scale)), interpolation=cv2.INTER_AREA)
    img_norm = normalize_image(img)
    cx, cy, r = find_ring_radius(img_norm)
    I_r_theta, radii = sample_I_r_theta(img_norm, cx, cy, r, r_fractions)

    # real
    eta_odd, _ = sf.parity_odd_contrast_ring(I_r_theta)
    eta_spiral = sf.spiral_harmonic_contrast(I_r_theta)
    knee_res = sf.knee_on_loglog(radii, eta_odd) if np.all(np.isfinite(eta_odd)) else None

    # null (theta shuffle)
    I_null = theta_shuffle(I_r_theta, rng)
    eta_odd_null, _ = sf.parity_odd_contrast_ring(I_null)
    eta_spiral_null = sf.spiral_harmonic_contrast(I_null)
    knee_null = sf.knee_on_loglog(radii, eta_odd_null) if np.all(np.isfinite(eta_odd_null)) else None

    return {
        "file": path.name,
        "r_base": r,
        "cx": cx,
        "cy": cy,
        "eta_odd_med": float(np.median(eta_odd)),
        "eta_odd_knee": knee_res.r_knee if knee_res else np.nan,
        "eta_spiral_med": float(np.median(eta_spiral)),
        "eta_odd_med_null": float(np.median(eta_odd_null)),
        "eta_odd_knee_null": knee_null.r_knee if knee_null else np.nan,
        "eta_spiral_med_null": float(np.median(eta_spiral_null)),
    }


def main(max_images: int = 40):
    rng = np.random.default_rng(42)
    roots = [
        Path("Additional data/converted_png"),
        Path("Additional data/converted_png_5761427"),
    ]
    r_fractions = np.linspace(0.6, 0.9, 4)
    rows = []
    for root in roots:
        if not root.exists():
            continue
        for path in sorted(root.glob("*.png")):
            if len(rows) >= max_images:
                break
            try:
                res = process_image(path, r_fractions, rng)
                res["dataset"] = root.name
                rows.append(res)
            except Exception as e:
                print(f"skip {path}: {e}")
                continue
    if not rows:
        raise SystemExit("No images processed")
    df = pd.DataFrame(rows)
    out_path = Path("Results/summaries/symmetry_meta.csv")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    print(f"Wrote {out_path}")
    print(df.head())


if __name__ == "__main__":
    main()
