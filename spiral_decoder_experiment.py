"""
Synthetic spiral/OAM encoding-decoding experiment with blind controls.

This script builds a challenge–response style benchmark:
1) Generate simple binary/greyscale patterns.
2) Encode them with a Fourier-domain spiral phase mask (or decoys).
3) Attempt recovery with a grid of candidate spiral parameters.
4) Compare to wrong-geometry decoders and generic enhancement baselines.
5) Log metrics (SSIM/PSNR/NCC) for matched vs control cases.

Usage (typical):
    python spiral_decoder_experiment.py --out-dir Results/spiral_trials --trials 8
    python spiral_decoder_experiment.py --out-dir Results/spiral_trials --no-save-images
"""
import argparse
import json
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple

import cv2
import numpy as np
import pandas as pd
from numpy.fft import fft2, ifft2


# --------------------------
# Pattern generation helpers
# --------------------------
def generate_patterns(size: int = 256) -> Dict[str, np.ndarray]:
    """
    Create a small library of reference patterns (float32 in [0,1]).
    """
    patterns: Dict[str, np.ndarray] = {}
    canvas = np.zeros((size, size), np.float32)

    def add_pattern(name: str, img: np.ndarray) -> None:
        patterns[name] = img.astype(np.float32)

    # Squares / bars
    square = canvas.copy()
    margin = size // 6
    cv2.rectangle(
        square,
        (margin, margin),
        (size - margin, size - margin),
        color=1.0,
        thickness=-1,
    )
    add_pattern("square", square)

    bars = canvas.copy()
    for k in range(0, size, size // 8):
        cv2.rectangle(bars, (k, 0), (k + size // 16, size), color=1.0, thickness=-1)
    add_pattern("bars", bars)

    checker = canvas.copy()
    block = size // 8
    for y in range(0, size, block):
        for x in range(0, size, block):
            if (x // block + y // block) % 2 == 0:
                cv2.rectangle(
                    checker, (x, y), (x + block, y + block), color=1.0, thickness=-1
                )
    add_pattern("checker", checker)

    ring = canvas.copy()
    cv2.circle(ring, (size // 2, size // 2), size // 3, color=1.0, thickness=8)
    add_pattern("ring", ring)

    spiral = canvas.copy()
    center = (size // 2, size // 2)
    angle = 0.0
    radius = size // 12
    while radius < size // 2:
        x = int(center[0] + radius * np.cos(angle))
        y = int(center[1] + radius * np.sin(angle))
        cv2.circle(spiral, (x, y), 3, color=1.0, thickness=-1)
        angle += np.pi / 10
        radius += 2
    spiral = cv2.GaussianBlur(spiral, (0, 0), 1.0)
    add_pattern("spiral_hint", spiral)

    letters = canvas.copy()
    cv2.putText(
        letters,
        "OAM",
        (size // 8, size // 2),
        cv2.FONT_HERSHEY_SIMPLEX,
        1.5,
        color=1.0,
        thickness=3,
        lineType=cv2.LINE_AA,
    )
    add_pattern("letters", letters)

    comma = canvas.copy()
    cv2.ellipse(comma, (size // 2, size // 2), (size // 4, size // 6), 20, 200, 340, color=1.0, thickness=5)
    cv2.circle(comma, (int(size * 0.65), int(size * 0.65)), size // 20, color=1.0, thickness=-1)
    add_pattern("comma", cv2.GaussianBlur(comma, (0, 0), 1.2))

    return patterns


# --------------------------
# Phase masks and encoding
# --------------------------
@dataclass
class MaskParams:
    kind: str  # "spiral", "random", "flat"
    ell: int = 1
    chirality: int = 1  # +1 or -1
    pitch: float = 0.0
    knee: float = 0.0
    seed: Optional[int] = None


@dataclass
class NoiseParams:
    gaussian_sigma: float = 0.01
    poisson_scale: float = 0.0


def make_phase_mask(shape: Tuple[int, int], params: MaskParams) -> np.ndarray:
    """Return complex phase mask exp(i*phi) with the requested geometry."""
    ny, nx = shape
    y = np.linspace(-1.0, 1.0, ny)
    x = np.linspace(-1.0, 1.0, nx)
    X, Y = np.meshgrid(x, y)
    R = np.sqrt(X ** 2 + Y ** 2) + 1e-6
    # Image coords have y increasing downward; flip Y so ell sign matches the helical projection convention.
    theta = np.arctan2(-Y, X)

    if params.kind == "spiral":
        # Spiral/OAM-informed phase: ell*theta plus optional radial terms.
        phi = params.chirality * params.ell * theta + params.pitch * R + params.knee * np.log(
            1.0 + R
        )
    elif params.kind == "random":
        rng = np.random.default_rng(params.seed)
        phi = rng.uniform(-np.pi, np.pi, size=shape)
    elif params.kind == "flat":
        phi = np.zeros_like(R)
    else:
        raise ValueError(f"Unknown mask kind: {params.kind}")

    return np.exp(1j * phi)


def normalize_image(img: np.ndarray) -> np.ndarray:
    img = np.clip(img, 0, None)
    mx = float(np.max(img))
    if mx == 0.0:
        return np.zeros_like(img, dtype=np.float32)
    return (img / mx).astype(np.float32)


def apply_noise(img: np.ndarray, noise: NoiseParams, rng: np.random.Generator) -> np.ndarray:
    noisy = img.copy()
    if noise.gaussian_sigma > 0:
        noisy = noisy + rng.normal(scale=noise.gaussian_sigma, size=img.shape).astype(
            np.float32
        )
    if noise.poisson_scale > 0:
        scaled = noisy * noise.poisson_scale
        noisy = rng.poisson(lam=np.clip(scaled, 0, None)).astype(np.float32)
        noisy = noisy / max(noise.poisson_scale, 1e-6)
    return normalize_image(noisy)


def encode_image(
    image: np.ndarray,
    mask_params: MaskParams,
    noise: NoiseParams,
    rng: np.random.Generator,
    ref_amp: float = 0.0,
    ref_kx: float = 0.0,
    ref_ky: float = 0.0,
) -> np.ndarray:
    """
    Encode an image by applying a phase mask in the Fourier domain and returning intensity.
    Optional: add an interferometric tilted reference wave to make phase observable.
    """
    # Apodize to limit spectral spread and reduce sideband overlap
    hann_y = np.hanning(image.shape[0])
    hann_x = np.hanning(image.shape[1])
    window = np.outer(hann_y, hann_x).astype(np.float32)
    image = image * window
    mask = make_phase_mask(image.shape, mask_params)
    F = fft2(image)
    F_masked = F * mask
    scrambled_field = ifft2(F_masked)
    if ref_amp != 0.0:
        h, w = image.shape
        y = np.linspace(-np.pi, np.pi, h)
        x = np.linspace(-np.pi, np.pi, w)
        X, Y = np.meshgrid(x, y)
        ref = ref_amp * np.exp(1j * (ref_kx * X + ref_ky * Y))
        scrambled_field = scrambled_field + ref
    intensity = np.abs(scrambled_field) ** 2
    intensity = normalize_image(intensity)
    return apply_noise(intensity, noise, rng)


def forward_intensity(
    image: np.ndarray, mask_params: MaskParams, ref_amp: float = 0.0, ref_kx: float = 0.0, ref_ky: float = 0.0
) -> np.ndarray:
    """
    Forward model without added noise (used for round-trip consistency scoring).
    """
    mask = make_phase_mask(image.shape, mask_params)
    F = fft2(image)
    F_masked = F * mask
    scrambled_field = ifft2(F_masked)
    if ref_amp != 0.0:
        h, w = image.shape
        y = np.linspace(-np.pi, np.pi, h)
        x = np.linspace(-np.pi, np.pi, w)
        X, Y = np.meshgrid(x, y)
        ref = ref_amp * np.exp(1j * (ref_kx * X + ref_ky * Y))
        scrambled_field = scrambled_field + ref
    intensity = np.abs(scrambled_field) ** 2
    return normalize_image(intensity)


# --------------------------
# Decoding and scoring
# --------------------------
def decode_with_mask(scrambled: np.ndarray, mask_params: MaskParams) -> np.ndarray:
    """
    Attempt to invert the mask: sqrt intensity -> FFT -> inverse phase -> IFFT.
    """
    mask = make_phase_mask(scrambled.shape, mask_params)
    amplitude = np.sqrt(np.clip(scrambled, 0, None))
    F = fft2(amplitude)
    restored = ifft2(F * np.conj(mask))
    decoded = np.real(restored)
    return normalize_image(decoded)


def detect_sideband(
    F_mag: np.ndarray,
    exclude_radius: int = 5,
    expand: float = 2.0,
) -> Tuple[Tuple[float, float], float, float, float]:
    """
    Detect strongest off-axis peak in |FFT|. Returns (y,x) center, radius, peak_snr, twin_leakage.
    """
    h, w = F_mag.shape
    cy, cx = h // 2, w // 2
    yy, xx = np.ogrid[:h, :w]
    center_mask = (yy - cy) ** 2 + (xx - cx) ** 2 <= exclude_radius ** 2
    masked = F_mag.copy()
    masked[center_mask] = 0.0
    idx = np.unravel_index(np.argmax(masked), masked.shape)
    peak_amp = masked[idx]
    win = 5
    y0, x0 = idx
    y_slice = slice(max(0, y0 - win), min(h, y0 + win + 1))
    x_slice = slice(max(0, x0 - win), min(w, x0 + win + 1))
    patch = masked[y_slice, x_slice]
    y_patch, x_patch = np.mgrid[y_slice, x_slice]
    weight = patch + 1e-9
    cy_sub = float(np.sum(y_patch * weight) / np.sum(weight))
    cx_sub = float(np.sum(x_patch * weight) / np.sum(weight))
    dy = y_patch - cy_sub
    dx = x_patch - cx_sub
    r2 = np.sum(weight * (dx ** 2 + dy ** 2)) / np.sum(weight)
    radius = max(exclude_radius + 2, expand * float(np.sqrt(max(r2, 1e-6))))
    median_amp = float(np.median(F_mag))
    peak_snr = float(peak_amp / (median_amp + 1e-9))
    twin_y = int(round(2 * cy - cy_sub))
    twin_x = int(round(2 * cx - cx_sub))
    twin_y = max(0, min(h - 1, twin_y))
    twin_x = max(0, min(w - 1, twin_x))
    twin_patch = F_mag[
        slice(max(0, twin_y - win), min(h, twin_y + win + 1)),
        slice(max(0, twin_x - win), min(w, twin_x + win + 1)),
    ]
    twin_leak = float(np.sum(twin_patch) / (np.sum(patch) + 1e-9))
    return (cy_sub, cx_sub), float(radius), peak_snr, twin_leak


def theta_image_coords(y: np.ndarray, x: np.ndarray, cy: float, cx: float) -> np.ndarray:
    """
    Angle with image coords (y down), so flip y to keep chirality consistent.
    """
    return np.arctan2(-(y - cy), (x - cx))


def local_peak(F_mag: np.ndarray, center: Tuple[float, float], window: int = 20) -> Tuple[float, float, float]:
    """
    Find peak in a local square window around center in fftshift coords.
    """
    h, w = F_mag.shape
    cy, cx = center
    y0 = int(max(0, cy - window))
    y1 = int(min(h, cy + window + 1))
    x0 = int(max(0, cx - window))
    x1 = int(min(w, cx + window + 1))
    patch = F_mag[y0:y1, x0:x1]
    if patch.size == 0:
        return cy, cx, 0.0
    idx = np.unravel_index(np.argmax(patch), patch.shape)
    py = y0 + idx[0]
    px = x0 + idx[1]
    val = float(patch[idx])
    return float(py), float(px), val


def adaptive_radius(
    F_mag: np.ndarray, center: Tuple[float, float], max_frac: float = 0.3, safety_margin: float = 4.0
) -> float:
    """
    Choose smallest radius capturing `max_frac` of sideband energy.
    """
    h, w = F_mag.shape
    cy, cx = center
    yy, xx = np.mgrid[0:h, 0:w]
    r = np.sqrt((yy - cy) ** 2 + (xx - cx) ** 2).flatten()
    energy = F_mag.flatten()
    order = np.argsort(r)
    r_sorted = r[order]
    e_sorted = energy[order]
    cumsum = np.cumsum(e_sorted)
    total = cumsum[-1] + 1e-9
    cutoff_idx = np.searchsorted(cumsum, max_frac * total)
    radius = float(r_sorted[min(cutoff_idx, len(r_sorted) - 1)] + safety_margin)
    return radius


def remove_phase_ramp(field: np.ndarray) -> np.ndarray:
    """
    Fit and remove a linear phase ramp (plane) from a complex field.
    """
    h, w = field.shape
    phase = np.angle(field)
    weight = np.abs(field)
    y, x = np.mgrid[0:h, 0:w]
    A = np.stack([x.flatten(), y.flatten(), np.ones(h * w)], axis=1)
    b = phase.flatten()
    w_flat = weight.flatten() + 1e-6
    Aw = A * w_flat[:, None]
    bw = b * w_flat
    coeffs, _, _, _ = np.linalg.lstsq(Aw, bw, rcond=None)
    plane = (coeffs[0] * x + coeffs[1] * y + coeffs[2])
    corrected = field * np.exp(-1j * plane.astype(np.float32))
    return corrected


def complex_log_polar(field: np.ndarray, center: Tuple[float, float], n_r: int, n_theta: int) -> np.ndarray:
    """
    Log-polar warp for complex field: warp real/imag separately then combine.
    """
    real = np.real(field).astype(np.float32)
    imag = np.imag(field).astype(np.float32)
    real_lp = cv2.warpPolar(real, (n_theta, n_r), center=center, maxRadius=np.hypot(center[0], center[1]), flags=cv2.WARP_POLAR_LOG + cv2.INTER_LINEAR + cv2.WARP_FILL_OUTLIERS)
    imag_lp = cv2.warpPolar(imag, (n_theta, n_r), center=center, maxRadius=np.hypot(center[0], center[1]), flags=cv2.WARP_POLAR_LOG + cv2.INTER_LINEAR + cv2.WARP_FILL_OUTLIERS)
    return real_lp.astype(np.float32) + 1j * imag_lp.astype(np.float32)


def oam_spectrum(field: np.ndarray, center: Tuple[float, float], n_r: int = 128, n_theta: int = 256, m_max: int = 10) -> Dict[str, float]:
    """
    Estimate OAM (azimuthal Fourier) spectrum power from complex field.
    Returns peak m and peakiness.
    """
    lp = complex_log_polar(field, center=center, n_r=n_r, n_theta=n_theta)
    # FFT over theta for each radius
    a_m = np.fft.fft(lp, axis=1)
    power = np.sum(np.abs(a_m) ** 2, axis=0)  # sum over radius
    power_shift = np.fft.fftshift(power)
    m_axis = np.arange(-n_theta // 2, n_theta // 2)
    # limit to |m|<=m_max neighborhood
    mask = np.abs(m_axis) <= m_max
    if not np.any(mask):
        return {"m_peak": 0.0, "m_peak_power": 0.0, "m_peak_ratio": 0.0}
    power_sel = power_shift[mask]
    m_sel = m_axis[mask]
    idx = int(np.argmax(power_sel))
    m_peak = float(m_sel[idx])
    peak_power = float(np.asarray(power_sel[idx]).item())
    total_power = float(np.sum(power_sel) + 1e-9)
    peak_ratio = float(peak_power / total_power)
    return {"m_peak": m_peak, "m_peak_power": peak_power, "m_peak_ratio": peak_ratio}


def carrier_from_ratio(radius: float, ratio: float) -> Tuple[float, float]:
    """
    Convert radius and direction ratio into (kx, ky) in cycles-per-image units.
    """
    ky = radius / max(np.sqrt(1.0 + ratio ** 2), 1e-6)
    kx = ratio * ky
    return kx, ky


def reconstruct_complex_field(
    hologram: np.ndarray,
    ref_kx: float,
    ref_ky: float,
    carrier_radius_frac: float = 0.12,
    auto_sideband: bool = True,
) -> Tuple[np.ndarray, Dict[str, float]]:
    """
    Off-axis holography demodulation: isolate sideband and shift to DC to recover complex field.
    If auto_sideband, detect peak location and radius; otherwise use ref_kx/ref_ky heuristics.
    """
    h, w = hologram.shape
    hann_y = np.hanning(h)
    hann_x = np.hanning(w)
    window = np.outer(hann_y, hann_x).astype(np.float32)
    hologram_win = hologram * window
    F = np.fft.fftshift(np.fft.fft2(hologram_win))
    cy, cx = h // 2, w // 2
    radius_init = max(5.0, carrier_radius_frac * min(h, w))
    r_min = 0.12 * min(h, w)
    min_crop_radius = max(6.0, 0.04 * min(h, w))
    fx_pred_raw = ref_kx / (2 * np.pi)
    fy_pred_raw = ref_ky / (2 * np.pi)
    fx_pred = ((fx_pred_raw + 0.5) % 1.0) - 0.5
    fy_pred = ((fy_pred_raw + 0.5) % 1.0) - 0.5
    peak_x_pred = cx + fx_pred * w
    peak_y_pred = cy + fy_pred * h
    r_max = 0.45 * min(h, w)
    kx_used_bins = fx_pred * w
    ky_used_bins = fy_pred * h
    def reconstruct_for_peak(peak_y: float, peak_x: float, radius: float) -> Tuple[np.ndarray, Dict[str, float]]:
        yy, xx = np.ogrid[:h, :w]
        twin_y = 2 * cy - peak_y
        twin_x = 2 * cx - peak_x
        max_allowed = np.sqrt((peak_y - cy) ** 2 + (peak_x - cx) ** 2) - 8.0
        radius_use = max(min_crop_radius, min(radius, max_allowed))
        mask = (yy - peak_y) ** 2 + (xx - peak_x) ** 2 <= radius_use ** 2
        mask_twin = (yy - twin_y) ** 2 + (xx - twin_x) ** 2 <= (radius_use ** 2)
        mask = np.logical_and(mask, np.logical_not(mask_twin))
        band = F * mask
        shift_x = int(round(cx - peak_x))
        shift_y = int(round(cy - peak_y))
        band_centered = np.roll(np.roll(band, shift_y, axis=0), shift_x, axis=1)
        field = np.fft.ifft2(np.fft.ifftshift(band_centered))
        patch = np.abs(F)[mask]
        twin_patch = np.abs(F)[mask_twin]
        # geometric overlap metric rather than energy parity of separate sidebands
        overlap = float(np.sum(mask & mask_twin) / (np.sum(mask) + 1e-9))
        dist_centers = float(np.hypot(peak_y - twin_y, peak_x - twin_x))
        leak_geom = max(0.0, 1.0 - dist_centers / (2.0 * radius_use + 1e-9))
        twin_leak = max(overlap, leak_geom)
        return field, {"crop_radius": float(radius_use), "twin_leak": float(twin_leak), "peak_y": float(peak_y), "peak_x": float(peak_x)}

    if auto_sideband:
        # local search around predicted lobes to avoid grabbing random object peaks
        py_plus, px_plus, _ = local_peak(np.abs(F), (peak_y_pred, peak_x_pred), window=24)
        py_minus, px_minus, _ = local_peak(np.abs(F), (2 * cy - peak_y_pred, 2 * cx - peak_x_pred), window=24)
        # also include a fallback global detection
        (peak_y_det, peak_x_det), radius_detected, peak_snr, _ = detect_sideband(np.abs(F))
        radius_energy = adaptive_radius(np.abs(F), (peak_y_det, peak_x_det), max_frac=0.8, safety_margin=4.0)
        peaks_to_eval = [
            (py_plus, px_plus),
            (py_minus, px_minus),
            (peak_y_det, peak_x_det),
            (2 * cy - peak_y_det, 2 * cx - peak_x_det),
        ]
        # evaluate candidate peaks and pick lowest leak
        candidates = []
        for py, px in peaks_to_eval:
            dist_dc = np.hypot(py - cy, px - cx)
            # gate: require proximity to predicted ±k
            gate = 20.0
            dist_plus = np.hypot(py - peak_y_pred, px - peak_x_pred)
            dist_minus = np.hypot(py - (2 * cy - peak_y_pred), px - (2 * cx - peak_x_pred))
            if min(dist_plus, dist_minus) > gate:
                continue
            if dist_dc < r_min or dist_dc > r_max:
                continue
            radius_energy_cand = adaptive_radius(np.abs(F), (py, px), max_frac=0.3, safety_margin=4.0)
            dist_twin = np.hypot(py - (2 * cy - py), px - (2 * cx - px))
            radius_geom_cap = 0.45 * max(1.0, min(dist_dc, dist_twin))
            radius_base = min(
                radius_geom_cap,
                r_max,
                max(radius_init, radius_detected, radius_energy_cand, min_crop_radius, dist_dc / 2.0 - 2.0),
            )
            # radius sweep to reduce leak
            radius_try = radius_base
            best_field = None
            best_diag = None
            for _ in range(6):
                field_cand, diag_cand = reconstruct_for_peak(py, px, radius_try)
                candidates.append((diag_cand["twin_leak"], diag_cand, field_cand))
                if diag_cand["twin_leak"] < 0.2:
                    break
                radius_try = max(min_crop_radius, radius_try * 0.8)
        if not candidates:
            # Fallback: use predicted peak directly
            field, diag = reconstruct_for_peak(peak_y_pred, peak_x_pred, radius_init)
            twin_leak = diag["twin_leak"]
        else:
            candidates.sort(key=lambda t: t[0])
            twin_leak = candidates[0][0]
            diag = candidates[0][1]
            field = candidates[0][2]
        # distances to predicted lobes for the chosen candidate
        dist_plus = float(np.hypot(diag["peak_y"] - peak_y_pred, diag["peak_x"] - peak_x_pred))
        dist_minus = float(np.hypot(diag["peak_y"] - (2 * cy - peak_y_pred), diag["peak_x"] - (2 * cx - peak_x_pred)))
        diag.update({"sideband_snr": float(peak_snr)})
        diag.update(
            {
                "peak_x_pred": float(peak_x_pred),
                "peak_y_pred": float(peak_y_pred),
                "r_min": float(r_min),
                "kx_used_bins": float(kx_used_bins),
                "ky_used_bins": float(ky_used_bins),
                "dist_detect_pred_plus": dist_plus,
                "dist_detect_pred_minus": dist_minus,
            }
        )
    else:
        peak_x = peak_x_pred
        peak_y = peak_y_pred
        radius_base = max(radius_init, min_crop_radius)
        field, diag = reconstruct_for_peak(peak_y, peak_x, radius_base)
        peak_snr = float(np.max(np.abs(F))) / (float(np.median(np.abs(F))) + 1e-9)
        diag.update({"sideband_snr": float(peak_snr)})
        twin_leak = diag["twin_leak"]
        diag.update(
            {
                "peak_x_pred": float(peak_x_pred),
                "peak_y_pred": float(peak_y_pred),
                "r_min": float(r_min),
                "kx_used_bins": float(kx_used_bins),
                "ky_used_bins": float(ky_used_bins),
                "dist_detect_pred_plus": float(np.hypot(peak_y - peak_y_pred, peak_x - peak_x_pred)),
                "dist_detect_pred_minus": float(np.hypot(peak_y - (2 * cy - peak_y_pred), peak_x - (2 * cx - peak_x_pred))),
            }
        )
    diag["twin_leak_raw"] = float(twin_leak)
    diag["twin_leak"] = float(np.clip(twin_leak, 0.0, 10.0))
    return field, diag


def decode_offaxis(
    hologram: np.ndarray,
    mask_params: MaskParams,
    ref_kx: float,
    ref_ky: float,
    carrier_radius_frac: float = 0.12,
    auto_sideband: bool = True,
) -> Tuple[np.ndarray, Dict[str, float]]:
    """
    Decode using off-axis hologram: demodulate to complex field, then invert mask.
    """
    field, diag = reconstruct_complex_field(
        hologram,
        ref_kx=ref_kx,
        ref_ky=ref_ky,
        carrier_radius_frac=carrier_radius_frac,
        auto_sideband=auto_sideband,
    )
    # If the selected lobe is the conjugate of the injected carrier, flip via conjugation to keep handedness consistent.
    dist_plus = diag.get("dist_detect_pred_plus", float("inf"))
    dist_minus = diag.get("dist_detect_pred_minus", float("inf"))
    if dist_minus < dist_plus:
        field = np.conj(field)
    mask = make_phase_mask(hologram.shape, mask_params)
    field = remove_phase_ramp(field)
    restored = np.fft.ifft2(np.fft.fft2(field) * np.conj(mask))
    decoded = np.real(restored)
    return normalize_image(decoded), diag


def phase_winding(field: np.ndarray, center: Tuple[float, float]) -> float:
    """
    Estimate phase winding number around a ring (independent chirality/ℓ witness).
    """
    h, w = field.shape
    radii = np.linspace(0.25, 0.45, 5) * (min(h, w) / 2.0)
    angles = np.linspace(-np.pi, np.pi, 256, endpoint=False)
    windings = []
    for r in radii:
        # Image coords: y increases downward, so flip sin to preserve handedness.
        xs = center[0] + r * np.cos(angles)
        ys = center[1] - r * np.sin(angles)
        xs = np.clip(xs, 0, w - 1)
        ys = np.clip(ys, 0, h - 1)
        phase_samples = []
        amp_samples = []
        for x, y in zip(xs, ys):
            i0, j0 = int(round(y)), int(round(x))
            val = field[i0, j0]
            phase_samples.append(np.angle(val))
            amp_samples.append(np.abs(val))
        amp_samples = np.array(amp_samples)
        if np.median(amp_samples) <= 1e-6:
            continue
        phase_unwrapped = np.unwrap(np.array(phase_samples))
        total = phase_unwrapped[-1] - phase_unwrapped[0]
        windings.append(total / (2 * np.pi))
    if not windings:
        return float("nan")
    return float(np.median(windings))


def estimate_vortex_center(
    field: np.ndarray, inner_margin_frac: float = 0.2, amp_weight: float = 0.5
) -> Tuple[float, float]:
    """
    Robust vortex core finder using phase residue (topological charge) with bounds.
    Returns (cx, cy) in image coordinates.
    """
    h, w = field.shape
    # Unit phasor
    z = field / (np.abs(field) + 1e-9)
    phi = np.angle(z)

    def wrap_pi(a: np.ndarray) -> np.ndarray:
        return (a + np.pi) % (2 * np.pi) - np.pi

    dpx = wrap_pi(phi[:, 1:] - phi[:, :-1])
    dpy = wrap_pi(phi[1:, :] - phi[:-1, :])
    circ = wrap_pi(dpx[:-1, :] + dpy[:, 1:] - dpx[1:, :] - dpy[:, :-1])
    charge = circ / (2 * np.pi)  # roughly integer per plaquette

    # Weight charge magnitude by smoothed amplitude to avoid noise domination
    amp = cv2.GaussianBlur(np.abs(field).astype(np.float32), (0, 0), 1.5)
    amp_pad = amp[1:, 1:]  # align with charge grid (h-1, w-1)
    score = np.abs(charge) * (amp_pad ** amp_weight)

    # Inner window to avoid edge picks
    y0 = int(inner_margin_frac * (h - 1))
    y1 = int((1 - inner_margin_frac) * (h - 1))
    x0 = int(inner_margin_frac * (w - 1))
    x1 = int((1 - inner_margin_frac) * (w - 1))
    if y1 <= y0 or x1 <= x0:
        y0, x0, y1, x1 = 0, 0, h - 1, w - 1
    mask = np.zeros_like(score, dtype=bool)
    mask[y0:y1, x0:x1] = True
    if np.any(mask):
        score_masked = np.where(mask, score, -np.inf)
    else:
        score_masked = score

    idx = int(np.argmax(score_masked))
    cy_cell, cx_cell = np.unravel_index(idx, score_masked.shape)
    # cell (cy_cell, cx_cell) corresponds to pixel corner; take center of cell
    cx = cx_cell + 0.5
    cy = cy_cell + 0.5

    # Fallback to image center if selection is near edge or nan
    edge_margin = min(cx, cy, w - 1 - cx, h - 1 - cy)
    if not np.isfinite(cx) or not np.isfinite(cy) or edge_margin < max(4.0, 0.05 * min(h, w)):
        cx = w / 2.0
        cy = h / 2.0
    return float(cx), float(cy)


def charge_metrics(
    field: np.ndarray,
    inner_margin_frac: float = 0.05,
    amp_percentile: float = 60.0,
    center: Optional[Tuple[float, float]] = None,
    r_inner_frac: float = 0.08,
    r_outer_frac: float = 0.4,
    eps: float = 1e-9,
) -> Dict[str, float]:
    """
    Compute charge map diagnostics: count vortices, sum charge, centroid of |charge|.
    Uses plaquette residue (wrapped phase circulation) without early rounding.
    """
    assert np.iscomplexobj(field), "charge_metrics expects complex field"
    h, w = field.shape
    cx = (w - 1) / 2 if center is None else float(center[0])
    cy = (h - 1) / 2 if center is None else float(center[1])
    yy, xx = np.mgrid[0:h, 0:w]
    rr = np.hypot(xx - cx, yy - cy)
    r_min = r_inner_frac * (min(h, w) / 2.0)
    r_max = r_outer_frac * (min(h, w) / 2.0)
    annulus = (rr >= r_min) & (rr <= r_max)
    amp_full = np.abs(field)
    if np.any(annulus):
        amp_thr = float(np.percentile(amp_full[annulus], amp_percentile))
    else:
        amp_thr = float(np.percentile(amp_full, amp_percentile))
    z = field / (amp_full + eps)
    phi = np.angle(z)
    phase_std = float(np.std(phi))
    imag_max = float(np.max(np.abs(np.imag(field))))

    def wrap_pi(a: np.ndarray) -> np.ndarray:
        return (a + np.pi) % (2 * np.pi) - np.pi

    # Plaquette residue via wrapped edge differences (robust, no final wrap)
    p00 = phi[:-1, :-1]
    p10 = phi[1:, :-1]
    p11 = phi[1:, 1:]
    p01 = phi[:-1, 1:]
    d01 = wrap_pi(p01 - p00)  # top
    d12 = wrap_pi(p11 - p01)  # right
    d23 = wrap_pi(p10 - p11)  # bottom (11->10)
    d30 = wrap_pi(p00 - p10)  # left  (10->00)
    loop = d01 + d12 + d23 + d30  # in [-4π, 4π]
    q_float = loop / (2 * np.pi)  # continuous charge per cell

    # Invalidate cells that touch very low amplitude (phase unreliable)
    amp_cells = np.minimum.reduce(
        [
            amp_full[:-1, :-1],
            amp_full[1:, :-1],
            amp_full[:-1, 1:],
            amp_full[1:, 1:],
        ]
    )
    annulus_cells = annulus[:-1, :-1] & annulus[1:, :-1] & annulus[:-1, 1:] & annulus[1:, 1:]
    valid_amp = (amp_cells > max(1e-8, amp_thr)) & annulus_cells
    q_float = np.where(valid_amp, q_float, 0.0)

    y0 = int(inner_margin_frac * (h - 1))
    y1 = int((1 - inner_margin_frac) * (h - 1))
    x0 = int(inner_margin_frac * (w - 1))
    x1 = int((1 - inner_margin_frac) * (w - 1))
    q_win = q_float[y0:y1, x0:x1] if y1 > y0 and x1 > x0 else q_float

    # Significant charge mask based on integer rounding of w_float
    w_int = np.rint(q_float).astype(np.int32)
    w_int_win = w_int[y0:y1, x0:x1] if y1 > y0 and x1 > x0 else w_int
    sig_mask = (np.abs(w_int_win) >= 1) & np.isfinite(w_int_win)
    charge_abs = np.abs(q_win)
    n_vort = int(np.sum(sig_mask))
    sum_charge = float(np.sum(w_int_win))
    sum_charge_abs = float(np.sum(np.abs(w_int_win)))
    q_mean = float(np.mean(q_win)) if q_win.size else float("nan")
    q_std = float(np.std(q_win)) if q_win.size else float("nan")
    q_max_abs = float(np.max(np.abs(q_win))) if q_win.size else float("nan")
    n_q_gt_0p25 = int(np.sum(np.abs(q_win) > 0.25))
    n_q_gt_0p5 = int(np.sum(np.abs(q_win) > 0.5))
    valid_win = valid_amp[y0:y1, x0:x1] if y1 > y0 and x1 > x0 else valid_amp
    q_valid = np.abs(q_win[valid_win]) if np.any(valid_win) else np.array([])
    q_abs_p50 = float(np.percentile(q_valid, 50)) if q_valid.size else float("nan")
    q_abs_p90 = float(np.percentile(q_valid, 90)) if q_valid.size else float("nan")
    n_valid_cells = int(np.sum(valid_win)) if valid_win.size else 0
    valid_cell_frac = float(np.mean(valid_win)) if valid_win.size else float("nan")

    # Centroid weighted by |q| on significant cells; fallback handled by caller
    weight = np.where(sig_mask, np.abs(q_win), 0.0)
    yy, xx = np.mgrid[0:q_win.shape[0], 0:q_win.shape[1]]
    denom = np.sum(weight)
    if denom > 0:
        cx = float(np.sum((xx + x0) * weight) / denom)
        cy = float(np.sum((yy + y0) * weight) / denom)
    else:
        cx = float("nan")
        cy = float("nan")

    return {
        "n_vort": float(n_vort),
        "sum_charge": sum_charge,
        "sum_charge_abs": sum_charge_abs,
        "q_float_mean": q_mean,
        "q_float_std": q_std,
        "q_float_max_abs": q_max_abs,
        "n_q_gt_0p25": float(n_q_gt_0p25),
        "n_q_gt_0p5": float(n_q_gt_0p5),
        "q_abs_p50": q_abs_p50,
        "q_abs_p90": q_abs_p90,
        "n_valid_cells": float(n_valid_cells),
        "valid_cell_frac": valid_cell_frac,
        "charge_cx": cx,
        "charge_cy": cy,
        "phase_std": phase_std,
        "imag_max": imag_max,
    }


def ensure_phase_alive(field: np.ndarray, tag: str = "") -> Tuple[float, float]:
    """
    Tripwire to catch accidentally real/flat fields before topology diagnostics.
    Raises RuntimeError if the phase looks dead.
    """
    if not np.iscomplexobj(field):
        raise RuntimeError(f"{tag}: field is not complex (dtype={getattr(field,'dtype',None)})")
    imag_max = float(np.max(np.abs(np.imag(field))))
    phase_std = float(np.std(np.angle(field)))
    if (not np.isfinite(phase_std)) or (imag_max < 1e-8) or (phase_std < 1e-4):
        raise RuntimeError(
            f"{tag}: phase looks dead (imag_max={imag_max:.3e}, phase_std={phase_std:.3e})"
        )
    return phase_std, imag_max


def field_fingerprint(field: np.ndarray) -> Dict[str, float]:
    """
    Lightweight fingerprint of a complex field to catch swaps/tap mistakes.
    """
    re = np.real(field)
    im = np.imag(field)
    return {
        "mean_real": float(np.mean(re)),
        "std_real": float(np.std(re)),
        "mean_imag": float(np.mean(im)),
        "std_imag": float(np.std(im)),
        "mean_abs": float(np.mean(np.abs(field))),
    }


def oam_spectrum_annulus(
    field: np.ndarray,
    center: Tuple[float, float],
    r_inner_frac: float = 0.08,
    r_outer_frac: float = 0.4,
    n_rings: int = 24,
    n_theta: int = 256,
    m_max: int = 3,
) -> Dict[str, float]:
    """
    OAM spectrum on an annulus using multiple rings; median across rings.
    """
    h, w = field.shape
    r_min = r_inner_frac * (min(h, w) / 2.0)
    r_max = r_outer_frac * (min(h, w) / 2.0)
    radii = np.linspace(r_min, r_max, n_rings)
    power_accum = None
    valid = 0
    for r in radii:
        angles = np.linspace(-np.pi, np.pi, n_theta, endpoint=False)
        xs = center[0] + r * np.cos(angles)
        ys = center[1] - r * np.sin(angles)
        xs = xs.astype(np.float32)
        ys = ys.astype(np.float32)
        real = cv2.remap(np.real(field).astype(np.float32), xs, ys, cv2.INTER_LINEAR, borderMode=cv2.BORDER_REFLECT)
        imag = cv2.remap(np.imag(field).astype(np.float32), xs, ys, cv2.INTER_LINEAR, borderMode=cv2.BORDER_REFLECT)
        ring = (real + 1j * imag).reshape(-1)
        # amplitude weighting on the ring to suppress empty/noisy sectors
        amp = np.abs(ring)
        if np.median(amp) <= 1e-9:
            continue
        wgt = amp / (np.sum(amp) + 1e-9)
        ring = ring * wgt
        if np.median(np.abs(ring)) <= 1e-6:
            continue
        a_m = np.fft.fft(ring)
        power = np.abs(a_m) ** 2
        if power_accum is None:
            power_accum = power
        else:
            power_accum += power
        valid += 1
    if valid == 0 or power_accum is None:
        return {
            "m_peak": 0.0,
            "m_peak_power": 0.0,
            "m_peak_ratio": 0.0,
            "total_power": 0.0,
        }
    power_accum = power_accum / valid
    power_shift = np.fft.fftshift(power_accum)
    m_axis = np.arange(-n_theta // 2, n_theta // 2)
    mask = np.abs(m_axis) <= m_max
    power_sel = power_shift[mask]
    m_sel = m_axis[mask]
    if power_sel.size == 0:
        return {
            "m_peak": 0.0,
            "m_peak_power": 0.0,
            "m_peak_ratio": 0.0,
            "total_power": 0.0,
        }
    idx = int(np.nanargmax(np.where(np.isfinite(power_sel), power_sel, -np.inf)))
    m_peak = float(m_sel[idx])
    peak_power = float(np.asarray(power_sel[idx]).item())
    total_power = float(np.sum(power_sel) + 1e-9)
    peak_ratio = float(peak_power / total_power if total_power > 0 else 0.0)
    return {
        "m_peak": m_peak,
        "m_peak_power": peak_power,
        "m_peak_ratio": peak_ratio,
        "total_power": float(total_power),
    }


def helical_projection(
    field: np.ndarray,
    center: Tuple[float, float],
    ells: Sequence[int] = (-3, -2, -1, 1, 2, 3),
    r_inner_frac: float = 0.08,
    r_outer_frac: float = 0.4,
    p: float = 1.0,
    n_rings: int = 1,
    phase_only: bool = False,
    scan_offsets: Sequence[Tuple[float, float]] = ((0.0, 0.0),),
    eps: float = 1e-12,
) -> Dict[str, float]:
    """
    Project the complex field onto helical basis exp(-i * ell * theta) over an annulus.
    Supports multi-annulus averaging, phase-only projection, and tiny center scans.
    Returns magnitudes per ell in [0,1] plus a few diagnostic fields (proj_peak_*).
    """
    h, w = field.shape
    yy, xx = np.mgrid[0:h, 0:w]
    r_min = r_inner_frac * (min(h, w) / 2.0)
    r_max = r_outer_frac * (min(h, w) / 2.0)
    if n_rings < 1:
        n_rings = 1
    r_edges = np.linspace(r_min, r_max, n_rings + 1)

    # Precompute amplitude; phase-only uses the unit phasor.
    amp_full = np.abs(field)
    if phase_only:
        z_full = field / (amp_full + eps)
    else:
        z_full = field

    def project_for_center(cx: float, cy: float) -> Tuple[Dict[str, float], float, float, float]:
        dx = xx - cx
        dy = cy - yy  # image coords (y down)
        r = np.hypot(dx, dy)
        th = theta_image_coords(yy, xx, cy, cx)
        num_acc = {ell: 0.0 for ell in ells}
        denom_total = 0.0
        for i in range(n_rings):
            r0, r1 = r_edges[i], r_edges[i + 1]
            ann = (r >= r0) & (r <= r1 if i == n_rings - 1 else r < r1)
            if not np.any(ann):
                continue
            z = z_full[ann]
            a = amp_full[ann]
            wgt = (a + eps) ** p
            denom = float(np.sum(a * wgt) + eps)
            denom_total += denom
            th_ann = th[ann]
            for ell in ells:
                basis = np.exp(-1j * ell * th_ann)
                num = np.abs(np.sum(z * basis * wgt))
                num_acc[ell] += float(num)
        if denom_total <= 0.0:
            return {f"proj_ell_{ell}": 0.0 for ell in ells}, 0.0, 0.0, 0.0
        out = {f"proj_ell_{ell}": num_acc[ell] / denom_total for ell in ells}
        vals = [(ell, out[f"proj_ell_{ell}"]) for ell in ells]
        vals.sort(key=lambda t: t[1], reverse=True)
        peak_ell, peak_val = vals[0]
        second_val = vals[1][1] if len(vals) > 1 else eps
        ratio = peak_val / (second_val + eps)
        return out, float(peak_ell), float(peak_val), float(ratio)

    best_proj: Dict[str, float] = {}
    best_ratio = -np.inf
    best_peak_val = -np.inf
    best_peak_ell = 0.0
    best_offset = (0.0, 0.0)
    for dx_off, dy_off in scan_offsets:
        proj, peak_ell, peak_val, ratio = project_for_center(center[0] + dx_off, center[1] + dy_off)
        # Prefer higher ratio; break ties with peak magnitude.
        if (ratio > best_ratio + 1e-9) or (np.isclose(ratio, best_ratio) and peak_val > best_peak_val):
            best_proj = proj
            best_ratio = ratio
            best_peak_val = peak_val
            best_peak_ell = peak_ell
            best_offset = (dx_off, dy_off)

    if not best_proj:
        best_proj = {f"proj_ell_{ell}": 0.0 for ell in ells}
        best_ratio = 0.0
        best_peak_val = 0.0
        best_peak_ell = 0.0
    best_proj.update(
        {
            "proj_peak_ell": float(best_peak_ell),
            "proj_peak_val": float(best_peak_val),
            "proj_second_val": float(best_peak_val / (best_ratio + eps)) if best_ratio > 0 else 0.0,
            "proj_peak_ratio": float(best_ratio),
            "proj_center_dx": float(best_offset[0]),
            "proj_center_dy": float(best_offset[1]),
            "proj_n_rings": int(n_rings),
            "proj_phase_only": int(bool(phase_only)),
        }
    )
    return best_proj


def helical_spectrum_stability(
    field: np.ndarray,
    center: Tuple[float, float],
    ells: Sequence[int] = (-3, -2, -1, 1, 2, 3),
    bands: Sequence[Tuple[float, float]] = ((0.05, 0.45), (0.05, 0.25), (0.25, 0.45)),
    center_offsets: Sequence[Tuple[float, float]] = ((0.0, 0.0), (1.0, 0.0), (-1.0, 0.0), (0.0, 1.0), (0.0, -1.0)),
) -> Dict[str, float]:
    """
    Robust spectrum estimate: aggregate phase-only helical projections across annuli and small center jitters.
    Returns median/std per ell plus peak/stability diagnostics.
    """
    records: List[Dict[str, float]] = []
    for r0, r1 in bands:
        for dx, dy in center_offsets:
            proj = helical_projection(
                field,
                center=(center[0] + dx, center[1] + dy),
                ells=ells,
                r_inner_frac=r0,
                r_outer_frac=r1,
                p=0.0,
                n_rings=3,
                phase_only=True,
                scan_offsets=((0.0, 0.0),),
            )
            records.append(proj)
    if not records:
        return {"spec_peak_ell": 0.0, "spec_peak_ratio": 0.0, "spec_stability_frac": 0.0, "spec_n": 0}
    vals_per_ell: Dict[int, List[float]] = {ell: [] for ell in ells}
    peak_ells: List[int] = []
    peak_ratios: List[float] = []
    for rec in records:
        for ell in ells:
            vals_per_ell[ell].append(float(rec.get(f"proj_ell_{ell}", 0.0)))
        peak_ells.append(int(rec.get("proj_peak_ell", 0)))
        peak_ratios.append(float(rec.get("proj_peak_ratio", 0.0)))
    med_spec = {ell: float(np.median(vals)) for ell, vals in vals_per_ell.items()}
    std_spec = {ell: float(np.std(vals)) for ell, vals in vals_per_ell.items()}
    med_items = sorted([(ell, val) for ell, val in med_spec.items()], key=lambda t: t[1], reverse=True)
    peak_ell = med_items[0][0]
    peak_val = med_items[0][1]
    second_val = med_items[1][1] if len(med_items) > 1 else 1e-9
    peak_ratio = float(peak_val / (second_val + 1e-9))
    stability = float(peak_ells.count(peak_ell) / len(peak_ells))
    out = {
        "spec_peak_ell": float(peak_ell),
        "spec_peak_val": float(peak_val),
        "spec_peak_ratio": float(peak_ratio),
        "spec_stability_frac": stability,
        "spec_n": len(records),
        "spec_peak_ratio_med": float(np.median(peak_ratios)),
    }
    for ell in ells:
        out[f"spec_med_ell_{ell}"] = med_spec[ell]
        out[f"spec_std_ell_{ell}"] = std_spec[ell]
    return out


def phase_circulation_annulus(
    field: np.ndarray,
    center: Tuple[float, float],
    r_inner_frac: float = 0.08,
    r_outer_frac: float = 0.4,
    n_rings: int = 24,
    n_theta: int = 256,
) -> Dict[str, float]:
    """
    Estimate topological charge via wrapped phase circulation over rings.
    """
    h, w = field.shape
    r_min = r_inner_frac * (min(h, w) / 2.0)
    r_max = r_outer_frac * (min(h, w) / 2.0)
    radii = np.linspace(r_min, r_max, n_rings)
    charges = []
    kept = 0
    total = 0
    for r in radii:
        angles = np.linspace(-np.pi, np.pi, n_theta, endpoint=False)
        xs = center[0] + r * np.cos(angles)
        ys = center[1] - r * np.sin(angles)
        real = cv2.remap(np.real(field).astype(np.float32), xs.astype(np.float32), ys.astype(np.float32), cv2.INTER_LINEAR, borderMode=cv2.BORDER_REFLECT)
        imag = cv2.remap(np.imag(field).astype(np.float32), xs.astype(np.float32), ys.astype(np.float32), cv2.INTER_LINEAR, borderMode=cv2.BORDER_REFLECT)
        ring = (real + 1j * imag).reshape(-1)
        amp = np.abs(ring)
        total += len(amp)
        if np.median(amp) <= 1e-6:
            continue
        thr = np.percentile(amp, 60.0)
        msk = amp > max(1e-8, 0.5 * thr)
        if msk.sum() < 0.2 * len(amp):
            continue
        kept += int(msk.sum())
        theta = np.linspace(-np.pi, np.pi, n_theta, endpoint=False)[msk]
        phi = np.angle(ring[msk])
        amp_w = amp[msk]
        w = (amp_w ** 2) / (np.sum(amp_w ** 2) + 1e-9)
        phi_unwrap = np.unwrap(phi)
        t0 = np.average(theta, weights=w)
        p0 = np.average(phi_unwrap, weights=w)
        num = np.sum(w * (theta - t0) * (phi_unwrap - p0))
        den = np.sum(w * (theta - t0) ** 2) + 1e-12
        ell_hat = num / den if den > 0 else 0.0  # in radians per rad => approx ell
        charges.append(float(ell_hat))
    if not charges:
        return {"m_circ_med": 0.0, "m_circ_std": 0.0, "m_circ_n": 0, "m_circ_kept_frac": 0.0}
    return {
        "m_circ_med": float(np.median(charges)),
        "m_circ_std": float(np.std(charges)),
        "m_circ_n": len(charges),
        "m_circ_kept_frac": float(kept / total) if total > 0 else 0.0,
    }



def image_score(decoded: np.ndarray) -> float:
    """
    Blind validity prior: prefer binary-ish images with sparse, coherent edges.
    - Entropy penalty on histogram (encourages values near 0/1)
    - Total variation penalty (sparser edges win)
    - Edge coherence bonus (aligned gradients)
    """
    # Histogram entropy for 6 bins in [0,1]
    hist, _ = np.histogram(decoded, bins=6, range=(0.0, 1.0), density=True)
    hist = hist + 1e-6
    entropy = -float(np.sum(hist * np.log(hist)))

    # Total variation (L1 of gradient magnitudes)
    gx = cv2.Sobel(decoded, cv2.CV_32F, 1, 0, ksize=3)
    gy = cv2.Sobel(decoded, cv2.CV_32F, 0, 1, ksize=3)
    grad_mag = np.sqrt(gx ** 2 + gy ** 2)
    tv = float(np.mean(grad_mag))

    # Edge coherence: penalize isotropic gradients, reward dominant orientation
    theta = np.arctan2(gy, gx)
    coherence = float(np.abs(np.mean(np.exp(1j * theta))))

    # Lower entropy and TV are better; higher coherence is better.
    return (-entropy) - 0.5 * tv + 2.0 * coherence


def grid_decode_joint(
    scrambles: Sequence[np.ndarray],
    candidates: Sequence[MaskParams],
    center: Tuple[float, float],
    n_r: int,
    n_theta: int,
    obs_bands: Optional[List[Dict[str, float]]] = None,
    ref_amp: float = 0.0,
    ref_kx: float = 0.0,
    ref_ky: float = 0.0,
    ang_weight: float = 0.1,
    use_m_distance: bool = False,
    obs_features: Optional[Dict[str, float]] = None,
    physics_weight: float = 0.0,
    template_cache: Optional[Dict[Tuple[float, ...], Dict[str, float]]] = None,
    probe: Optional[np.ndarray] = None,
    decode_fn: Optional[Callable[[np.ndarray, MaskParams], np.ndarray]] = None,
    leak_penalty: float = 0.0,
    target_winding: Optional[float] = None,
    oam_target: Optional[Dict[str, float]] = None,
    oam_weight: float = 0.0,
) -> Tuple[np.ndarray, MaskParams, float, float, List[Dict[str, object]], float]:
    """
    Decode multiple scrambles jointly: score sums across inputs.
    Round-trip score: decode -> re-encode -> compare to observed (negative MSE).
    Adds angular feature consistency distance when observed bands are provided.
    Returns best decoded image for the first scramble, the best mask, best score, margin to runner-up,
    and top-5 scored candidates for diagnostics.
    """
    best_score = -np.inf
    second_score = -np.inf
    second_mask: Optional[MaskParams] = None
    best_mask = None
    best_img = None
    scored: List[Tuple[float, MaskParams, float, float, float]] = []
    best_physics = 0.0
    tmpl_cache = {} if template_cache is None else template_cache
    probe_img = probe
    if probe_img is None:
        probe_img = template_probe(scrambles[0].shape)
    decode_fn_use = decode_fn if decode_fn is not None else decode_with_mask

    def mask_equiv_key(p: MaskParams) -> Tuple[object, ...]:
        return (p.kind, p.chirality * p.ell, p.pitch, p.knee)

    def js_divergence(hist_a: Dict[int, int], hist_b: Dict[int, int]) -> float:
        keys = set(hist_a.keys()) | set(hist_b.keys())
        if not keys:
            return 0.0
        pa = np.array([hist_a.get(k, 0) for k in keys], dtype=np.float64)
        pb = np.array([hist_b.get(k, 0) for k in keys], dtype=np.float64)
        pa = pa / (pa.sum() + 1e-12)
        pb = pb / (pb.sum() + 1e-12)
        m = 0.5 * (pa + pb)
        def kl(p, q):
            mask = p > 0
            return float(np.sum(p[mask] * np.log((p[mask] + 1e-12) / (q[mask] + 1e-12))))
        return 0.5 * kl(pa, m) + 0.5 * kl(pb, m)

    def band_distance(obs: Dict[str, float], rec: Dict[str, float]) -> float:
        if obs.get("band_dead", 0) or rec.get("band_dead", 0):
            return 1.0
        d_fv = abs(obs.get("fv_peak_marg", 0.0) - rec.get("fv_peak_marg", 0.0))
        d_frac = abs(obs.get("energy_out", 0.0) / (obs.get("energy_total", 1e-9)) - rec.get("energy_out", 0.0) / (rec.get("energy_total", 1e-9)))
        d_ratio = abs(obs.get("fv_peak_ratio", 0.0) - rec.get("fv_peak_ratio", 0.0))
        d_js = js_divergence(obs.get("fv_hist_topk", {}), rec.get("fv_hist_topk", {}))
        d = d_fv + d_frac + d_ratio + d_js
        if use_m_distance:
            d_m = abs(obs.get("m_peak", 0.0) - rec.get("m_peak", 0.0))
            d_m_ratio = abs(obs.get("m_ratio", 0.0) - rec.get("m_ratio", 0.0))
            d += d_m + d_m_ratio
        return float(d)
    for params in candidates:
        total = 0.0
        decoded_first = None
        ang_dist_accum = []
        for idx, scrambled in enumerate(scrambles):
            decoded_out = decode_fn_use(scrambled, params)
            if isinstance(decoded_out, tuple):
                decoded = decoded_out[0]
            else:
                decoded = decoded_out
            if idx == 0:
                decoded_first = decoded
            reenc = forward_intensity(decoded, params, ref_amp=ref_amp, ref_kx=ref_kx, ref_ky=ref_ky)
            mse = float(np.mean((scrambled - reenc) ** 2))
            total += -mse
            if obs_bands is not None and len(scrambles) > 0:
                bands_rec = annular_band_peaks(log_polar_warp(reenc, center=center, n_r=n_r, n_theta=n_theta), n_bands=len(obs_bands))
                d_bands = []
                for b_obs, b_rec in zip(obs_bands, bands_rec):
                    d_bands.append(band_distance(b_obs, b_rec))
                if d_bands:
                    ang_dist_accum.append(float(np.median(d_bands)))

        ang_penalty = float(np.mean(ang_dist_accum)) if ang_dist_accum else 0.0
        physics_term = 0.0
        if physics_weight != 0.0 and obs_features is not None:
            key = (
                params.kind,
                params.ell,
                params.chirality,
                float(params.pitch),
                float(params.knee),
                float(n_r),
                float(n_theta),
            )
            if key not in tmpl_cache:
                tmpl_intensity = forward_intensity(probe_img, params, ref_amp=ref_amp, ref_kx=ref_kx, ref_ky=ref_ky)
                tmpl_cache[key] = spiral_logpolar_features(tmpl_intensity, center=center, n_r=n_r, n_theta=n_theta)
            tmpl_feat = tmpl_cache[key]
            physics_term = spiral_physics_likelihood(obs_features, tmpl_feat, params, n_theta=n_theta)

        # Theory-inspired gentle bias: prefer left-handed and larger |ell| only on ties
        tie_bias = 1e-8
        winding_pen = 0.0
        if target_winding is not None and np.isfinite(target_winding):
            winding_pen = abs(abs(target_winding) - abs(params.ell)) + (0.5 if np.sign(target_winding) != np.sign(params.chirality * params.ell) else 0.0)
        oam_term = 0.0
        if oam_target is not None:
            gate_reason = oam_target.get("gate_reason", "")
            proj_ratio = float(oam_target.get("proj_peak_ratio", oam_target.get("m_peak_ratio", 0.0)))
            proj_ratio_thresh = 1.1
            proj_weight = 0.05
            if gate_reason == "proj" and proj_ratio >= proj_ratio_thresh:
                reliability = np.clip(
                    (proj_ratio - proj_ratio_thresh) / max(proj_ratio_thresh, 1e-6),
                    0.0,
                    1.0,
                )
                c_pos = float(oam_target.get(f"proj_ell_{params.ell}", 0.0))
                c_neg = float(oam_target.get(f"proj_ell_{-params.ell}", 0.0))
                proj_diff = c_pos - c_neg
                oam_term = reliability * proj_weight * proj_diff
        total_biased = (
            total
            - ang_weight * ang_penalty
            + physics_weight * physics_term
            + oam_weight * oam_term
            - 1000.0 * leak_penalty
            - 0.5 * winding_pen
            + (tie_bias if params.chirality == -1 else 0.0)
            + 1e-6 * abs(params.ell)
        )
        scored.append((total_biased, params, physics_term, total, ang_penalty))
        if total_biased > best_score + 1e-12 or best_mask is None:
            second_score = best_score
            second_mask = best_mask
            best_score = total_biased
            best_mask = params
            best_img = decoded_first
            best_physics = physics_term
        elif total_biased > second_score + 1e-12:
            second_score = total_biased
            second_mask = params
    assert best_mask is not None and best_img is not None
    scored_sorted = sorted(scored, key=lambda x: x[0], reverse=True)
    best_entry = scored_sorted[0]
    runner_up_score = next(
        (s for s, p, _, _, _ in scored_sorted if mask_equiv_key(p) != mask_equiv_key(best_mask)),
        None,
    )
    margin = best_score - runner_up_score if runner_up_score is not None else float("nan")
    top5 = [
        {"score": float(s), "physics_score": float(p_physics), "roundtrip_score": float(rt), "ang_penalty": float(ap), **asdict(p)}
        for s, p, p_physics, rt, ap in scored_sorted[:5]
    ]
    return best_img, best_mask, float(best_score), float(margin), top5, float(best_physics)


# --------------------------
# Metrics
# --------------------------
def psnr(ref: np.ndarray, test: np.ndarray, eps: float = 1e-9) -> float:
    mse = float(np.mean((ref - test) ** 2))
    if mse == 0:
        return float("inf")
    return 10 * np.log10(1.0 / (mse + eps))


def ssim_simple(
    ref: np.ndarray, test: np.ndarray, c1: float = 0.01 ** 2, c2: float = 0.03 ** 2
) -> float:
    """
    Lightweight SSIM implementation (single-scale, grayscale).
    """
    mu_x = cv2.GaussianBlur(ref, (7, 7), 1.5)
    mu_y = cv2.GaussianBlur(test, (7, 7), 1.5)
    sigma_x = cv2.GaussianBlur(ref * ref, (7, 7), 1.5) - mu_x ** 2
    sigma_y = cv2.GaussianBlur(test * test, (7, 7), 1.5) - mu_y ** 2
    sigma_xy = cv2.GaussianBlur(ref * test, (7, 7), 1.5) - mu_x * mu_y
    numerator = (2 * mu_x * mu_y + c1) * (2 * sigma_xy + c2)
    denominator = (mu_x ** 2 + mu_y ** 2 + c1) * (sigma_x + sigma_y + c2)
    ssim_map = numerator / (denominator + 1e-9)
    return float(np.clip(np.mean(ssim_map), -1.0, 1.0))


def ncc(ref: np.ndarray, test: np.ndarray) -> float:
    ref_c = ref - float(np.mean(ref))
    test_c = test - float(np.mean(test))
    denom = np.sqrt(np.sum(ref_c ** 2) * np.sum(test_c ** 2)) + 1e-9
    return float(np.sum(ref_c * test_c) / denom)


# --------------------------
# Log-polar / annular analysis
# --------------------------
def log_polar_warp(image: np.ndarray, center: Tuple[float, float], n_r: int, n_theta: int) -> np.ndarray:
    """
    Map image to log-polar coordinates (u = log r, v = theta).
    """
    h, w = image.shape
    max_radius = np.hypot(max(center[0], w - center[0]), max(center[1], h - center[1]))
    log_base = np.exp(np.log(max_radius) / (n_r - 1))
    dst = cv2.warpPolar(
        image,
        (n_theta, n_r),
        center=center,
        maxRadius=max_radius,
        flags=cv2.WARP_POLAR_LOG + cv2.INTER_LINEAR + cv2.WARP_FILL_OUTLIERS,
    )
    return dst.astype(np.float32)


def estimate_center(img: np.ndarray) -> Tuple[float, float]:
    """
    Estimate image center via intensity centroid; fallback to geometric center.
    """
    h, w = img.shape
    total = float(np.sum(img))
    if total <= 0:
        return (w / 2.0, h / 2.0)
    y_idx, x_idx = np.mgrid[0:h, 0:w]
    cx = float(np.sum(x_idx * img) / total)
    cy = float(np.sum(y_idx * img) / total)
    return (cx, cy)


def annular_band_peaks(logpolar: np.ndarray, n_bands: int = 6, top_k: int = 50) -> List[Dict[str, float]]:
    """
    Split log-polar image into radial bands, FFT each, and extract dominant peak.
    Returns list of dicts with (fu, fv, amp).
    DC/near-DC bins are suppressed so angular structure (fv) can emerge.
    """
    bands = []
    n_r, n_theta = logpolar.shape
    edges = np.linspace(0, n_r, n_bands + 1, dtype=int)
    for i in range(n_bands):
        r0, r1 = edges[i], edges[i + 1]
        band = logpolar[r0:r1]
        if band.size == 0:
            bands.append({"fu": 0.0, "fv": 0.0, "amp": 0.0})
            continue
        # Zero-mean per-theta and detrend theta to suppress fv=0 ridge
        theta_mean = np.mean(band, axis=1, keepdims=True)
        theta_idx = np.arange(band.shape[1], dtype=np.float32)
        theta_idx = theta_idx - theta_idx.mean()
        theta_var = float(np.mean(np.var(band, axis=1)))  # pre-processed theta variance
        # Linear detrend per radius
        slope = (band @ theta_idx) / (np.sum(theta_idx ** 2) + 1e-9)
        band = band - theta_mean - slope[:, None] * theta_idx
        F = np.abs(np.fft.fftshift(np.fft.fft2(band)))
        # Suppress DC and near-DC neighborhood + wider fv corridor
        cx, cy = F.shape[0] // 2, F.shape[1] // 2
        win = 5
        F[cx - win : cx + win + 1, cy - win : cy + win + 1] = 0.0
        fv_half_width = max(3, int(0.05 * F.shape[1]))
        fv_mask = np.ones_like(F, dtype=bool)
        fv_mask[:, cy - fv_half_width : cy + fv_half_width + 1] = False
        masked = F * fv_mask
        energy_total = float(np.sum(F))
        energy_out = float(np.sum(masked))
        max_out = float(np.max(masked)) if masked.size else 0.0
        if max_out <= 1e-9 or energy_out <= 1e-9:
            bands.append(
                {
                    "fu": 0.0,
                    "fv": float("nan"),
                    "amp": 0.0,
                    "top_nonzero_fv": 0,
                    "top_abs_fv": 0.0,
                    "fv_mask_width": fv_half_width,
                    "band_dead": 1,
                    "energy_total": energy_total,
                    "energy_out": energy_out,
                    "fv_peak_marg": float("nan"),
                    "fv_peak_ratio": 0.0,
                    "fv_entropy": float("nan"),
                    "fv_hist_topk": [],
                    "theta_var": theta_var,
                    "m_peak": float("nan"),
                    "m_ratio": 0.0,
                    "m_entropy": float("nan"),
                }
            )
            continue

        flat_idx = np.argpartition(masked.flatten(), -top_k)[-top_k:]
        flat_sorted = flat_idx[np.argsort(masked.flatten()[flat_idx])[::-1]]

        def idx_to_freq(idx_val: int) -> Tuple[float, float]:
            i, j = np.unravel_index(idx_val, F.shape)
            fu = (i - F.shape[0] / 2) / max(F.shape[0], 1)
            fv = (j - F.shape[1] / 2) / max(F.shape[1], 1)
            return float(fu), float(fv)

        top_entries = []
        fv_hist = {}
        for k_idx in flat_sorted:
            amp = float(masked.flatten()[k_idx])
            fu, fv = idx_to_freq(int(k_idx))
            top_entries.append({"fu": fu, "fv": fv, "amp": amp})
            fv_bin = int(round(fv * F.shape[1]))
            fv_hist[fv_bin] = fv_hist.get(fv_bin, 0) + 1

        # Best peak with |fv|>0 due to mask
        best = top_entries[0]
        nonzero_fv = [e for e in top_entries if abs(e["fv"]) > 0]
        top_nonzero_fv = len(nonzero_fv)
        top_abs_fv = float(max(abs(e["fv"]) for e in top_entries)) if top_entries else 0.0

        # 1D marginal over fv (summed over fu)
        P_fv = masked.sum(axis=0)
        P_fv_sum = float(np.sum(P_fv))
        if P_fv_sum > 0:
            fv_star_idx = int(np.argmax(P_fv))
            fv_star = (fv_star_idx - F.shape[1] / 2) / max(F.shape[1], 1)
            fv_peak_ratio = float(P_fv[fv_star_idx] / (P_fv_sum + 1e-9))
            p_norm = P_fv / (P_fv_sum + 1e-9)
            fv_entropy = float(-np.sum(p_norm * np.log(p_norm + 1e-12)))
        else:
            fv_star = float("nan")
            fv_peak_ratio = 0.0
            fv_entropy = float("nan")

        # Circular harmonic over theta (collapse radius then FFT over theta)
        theta_signal = np.mean(band, axis=0)
        m_fft = np.abs(np.fft.fftshift(np.fft.fft(theta_signal)))
        m_half_width = max(1, int(0.05 * m_fft.size))
        m_fft[m_fft.size // 2 - m_half_width : m_fft.size // 2 + m_half_width + 1] = 0.0
        if np.sum(m_fft) > 0:
            m_star_idx = int(np.argmax(m_fft))
            m_star = (m_star_idx - m_fft.size / 2) / max(m_fft.size, 1)
            m_ratio = float(m_fft[m_star_idx] / (np.sum(m_fft) + 1e-9))
            m_norm = m_fft / (np.sum(m_fft) + 1e-9)
            m_entropy = float(-np.sum(m_norm * np.log(m_norm + 1e-12)))
        else:
            m_star = float("nan")
            m_ratio = 0.0
            m_entropy = float("nan")

        bands.append(
            {
                "fu": best["fu"],
                "fv": best["fv"],
                "amp": best["amp"],
                "top_nonzero_fv": top_nonzero_fv,
                "top_abs_fv": top_abs_fv,
                "fv_mask_width": fv_half_width,
                "fv_hist_topk": fv_hist,
                "band_dead": 0,
                "energy_total": energy_total,
                "energy_out": energy_out,
                "fv_peak_marg": fv_star,
                "fv_peak_ratio": fv_peak_ratio,
                "fv_entropy": fv_entropy,
                "theta_var": theta_var,
                "m_peak": m_star,
                "m_ratio": m_ratio,
                "m_entropy": m_entropy,
            }
        )
    return bands


def summarize_bands(bands: List[Dict[str, float]]) -> Dict[str, float]:
    """
    Turn per-band peaks into coarse G/S/E features.
    """
    if not bands:
        return {"fv_mean": 0.0, "fv_std": 0.0, "fu_mean": 0.0, "fu_std": 0.0, "amp_mean": 0.0, "amp_decay": 0.0}
    fv = np.array([b["fv"] for b in bands], dtype=np.float32)
    fu = np.array([b["fu"] for b in bands], dtype=np.float32)
    amp = np.array([b["amp"] for b in bands], dtype=np.float32)
    top_abs_fv = np.array([b.get("top_abs_fv", 0.0) for b in bands], dtype=np.float32)
    top_nonzero_fv = np.array([b.get("top_nonzero_fv", 0) for b in bands], dtype=np.float32)
    theta_var = np.array([b.get("theta_var", 0.0) for b in bands], dtype=np.float32)
    fv_peak_marg = np.array([b.get("fv_peak_marg", 0.0) for b in bands], dtype=np.float32)
    fv_peak_ratio = np.array([b.get("fv_peak_ratio", 0.0) for b in bands], dtype=np.float32)
    fv_entropy = np.array([b.get("fv_entropy", 0.0) for b in bands], dtype=np.float32)
    band_dead = np.array([b.get("band_dead", 0) for b in bands], dtype=np.float32)
    energy_out = np.array([b.get("energy_out", 0.0) for b in bands], dtype=np.float32)
    energy_total = np.array([b.get("energy_total", 0.0) for b in bands], dtype=np.float32)
    fv_nonzero_energy_frac = energy_out / (energy_total + 1e-9)
    fv_nonzero_energy_frac = np.clip(fv_nonzero_energy_frac, 0.0, 1.0)
    m_peak = np.array([b.get("m_peak", 0.0) for b in bands], dtype=np.float32)
    m_ratio = np.array([b.get("m_ratio", 0.0) for b in bands], dtype=np.float32)
    m_entropy = np.array([b.get("m_entropy", 0.0) for b in bands], dtype=np.float32)
    amp_decay = 0.0
    if len(amp) > 1:
        amp_decay = float(np.maximum(0.0, amp[0] - amp[-1])) / (float(amp[0]) + 1e-6)
    return {
        "fv_mean": float(np.mean(fv)),
        "fv_std": float(np.std(fv)),
        "fv_peak_marg_mean": float(np.nanmean(fv_peak_marg)),
        "fv_peak_marg_std": float(np.nanstd(fv_peak_marg)),
        "fv_peak_ratio_mean": float(np.nanmean(fv_peak_ratio)),
        "fv_entropy_mean": float(np.nanmean(fv_entropy)),
        "fv_nonzero_energy_frac_mean": float(np.nanmean(fv_nonzero_energy_frac)),
        "m_peak_mean": float(np.nanmean(m_peak)),
        "m_peak_std": float(np.nanstd(m_peak)),
        "m_ratio_mean": float(np.nanmean(m_ratio)),
        "m_entropy_mean": float(np.nanmean(m_entropy)),
        "fu_mean": float(np.mean(fu)),
        "fu_std": float(np.std(fu)),
        "amp_mean": float(np.mean(amp)),
        "amp_decay": amp_decay,
        "top_abs_fv_mean": float(np.mean(top_abs_fv)),
        "top_nonzero_fv_mean": float(np.mean(top_nonzero_fv)),
        "fv_mask_width_mean": float(np.mean([b.get("fv_mask_width", 0) for b in bands])),
        "theta_var_mean": float(np.mean(theta_var)),
        "band_dead_frac": float(np.mean(band_dead)),
    }


def angular_mode_features(logpolar: np.ndarray) -> Dict[str, float]:
    """
    Aggregate angular Fourier statistics from a log-polar image.
    """
    n_r, n_theta = logpolar.shape
    lp = logpolar - float(np.mean(logpolar))
    lp = cv2.GaussianBlur(lp, (3, 3), 0)
    m_fft = np.abs(np.fft.fftshift(np.fft.fft(lp, axis=1), axes=1))
    center = m_fft.shape[1] // 2
    dc = max(1, int(0.02 * m_fft.shape[1]))
    m_fft[:, center - dc : center + dc + 1] = 0.0
    spectrum = np.mean(m_fft, axis=0)
    total = float(np.sum(spectrum))
    if total <= 1e-9:
        return {
            "m_peak": 0.0,
            "m_ratio": 0.0,
            "m_entropy": float("nan"),
            "m_abs_mean": 0.0,
            "mode_energy": 0.0,
            "n_theta": n_theta,
        }
    spectrum = spectrum / (total + 1e-9)
    m_axis = (np.arange(len(spectrum), dtype=np.float32) - len(spectrum) / 2) / max(len(spectrum), 1)
    peak_idx = int(np.argmax(spectrum))
    m_peak = float(m_axis[peak_idx])
    m_ratio = float(np.max(spectrum))
    m_entropy = float(-np.sum(spectrum * np.log(spectrum + 1e-12)))
    m_abs_mean = float(np.sum(np.abs(m_axis) * spectrum))
    return {
        "m_peak": m_peak,
        "m_ratio": m_ratio,
        "m_entropy": m_entropy,
        "m_abs_mean": m_abs_mean,
        "mode_energy": total / float(n_r),
        "n_theta": n_theta,
    }


def ridge_slope_features(logpolar: np.ndarray) -> Dict[str, float]:
    """
    Estimate dominant ridge slope in log-polar space (pitch proxy).
    """
    blurred = cv2.GaussianBlur(logpolar, (5, 5), 0)
    gx = cv2.Sobel(blurred, cv2.CV_32F, 1, 0, ksize=3)
    gy = cv2.Sobel(blurred, cv2.CV_32F, 0, 1, ksize=3)
    mag = np.sqrt(gx ** 2 + gy ** 2)
    if not np.any(np.isfinite(mag)):
        return {"slope_med": 0.0, "slope_mean": 0.0, "slope_std": 0.0, "ridge_energy": 0.0}
    thresh = max(np.percentile(mag, 75), 1e-6)
    mask = mag > thresh
    ridge_orientation = np.arctan2(gy, gx) + np.pi / 2.0
    slope = np.tan(ridge_orientation)
    slope = np.clip(slope, -5.0, 5.0)
    slope_vals = slope[mask]
    mag_vals = mag[mask]
    if slope_vals.size == 0:
        return {"slope_med": 0.0, "slope_mean": 0.0, "slope_std": 0.0, "ridge_energy": 0.0}
    slope_med = float(np.median(slope_vals))
    slope_mean = float(np.sum(slope_vals * mag_vals) / (np.sum(mag_vals) + 1e-9))
    slope_std = float(np.std(slope_vals))
    ridge_energy = float(np.mean(mag_vals))
    return {"slope_med": slope_med, "slope_mean": slope_mean, "slope_std": slope_std, "ridge_energy": ridge_energy}


def change_point(signal: np.ndarray) -> Tuple[float, float]:
    """
    Lightweight change-point heuristic based on largest discrete derivative.
    """
    y = np.asarray(signal, dtype=np.float32).ravel()
    n = y.size
    if n < 4:
        return 0.0, float("inf")
    # Smooth to avoid noise-driven peaks
    kernel = np.ones(5, dtype=np.float32) / 5.0
    y_smooth = np.convolve(y, kernel, mode="same")
    dy = np.abs(np.diff(y_smooth))
    idx = int(np.argmax(dy))
    split_norm = idx / max(n - 1, 1)
    quality = float(np.mean(dy))
    return split_norm, quality


def knee_features(logpolar: np.ndarray) -> Dict[str, float]:
    """
    Detect knee/change-point in radial energy profile of angular content.
    """
    centered = logpolar - float(np.mean(logpolar))
    theta_fft = np.abs(np.fft.fftshift(np.fft.fft(centered, axis=1), axes=1))
    radial_energy = np.mean(theta_fft, axis=1)
    knee_idx, knee_err = change_point(radial_energy)
    return {
        "knee_idx": knee_idx,
        "knee_err": knee_err,
        "radial_energy_mean": float(np.mean(radial_energy)),
        "radial_energy_std": float(np.std(radial_energy)),
    }


def spiral_logpolar_features(image: np.ndarray, center: Tuple[float, float], n_r: int, n_theta: int) -> Dict[str, float]:
    """
    Convenience wrapper to compute angular, slope, and knee descriptors in log-polar space.
    """
    lp = log_polar_warp(image, center=center, n_r=n_r, n_theta=n_theta)
    ang = angular_mode_features(lp)
    ridge = ridge_slope_features(lp)
    knee = knee_features(lp)
    return {**ang, **ridge, **knee}


def template_probe(shape: Tuple[int, int]) -> np.ndarray:
    """
    Construct a simple probe pattern used for template feature prediction.
    """
    h, w = shape
    img = np.zeros((h, w), np.float32)
    c = (w // 2, h // 2)
    r = min(h, w) // 3
    cv2.circle(img, c, r, color=1.0, thickness=4)
    cv2.line(img, (c[0], c[1] - r), (c[0], c[1] + r), color=0.6, thickness=2)
    cv2.line(img, (c[0] - r, c[1]), (c[0] + r, c[1]), color=0.6, thickness=2)
    return normalize_image(cv2.GaussianBlur(img, (0, 0), 1.2))


def spiral_physics_likelihood(
    obs_feat: Dict[str, float],
    tmpl_feat: Dict[str, float],
    params: MaskParams,
    n_theta: int,
    m_tol: float = 0.02,
    slope_tol: float = 0.6,
    knee_tol: float = 0.15,
) -> float:
    """
    Compare observed log-polar features against a mask-informed template.
    Higher is better.
    """
    obs = {k: float(v) if np.isfinite(v) else 0.0 for k, v in obs_feat.items()}
    tmpl = {k: float(v) if np.isfinite(v) else 0.0 for k, v in tmpl_feat.items()}
    m_target = params.chirality * params.ell / max(float(n_theta), 1.0)
    m_obs = obs.get("m_peak", 0.0)
    m_tmp = tmpl.get("m_peak", 0.0)
    m_diff = abs(m_obs - m_target)
    m_consistency = abs(m_obs - m_tmp)
    mode_term = -(m_diff / (m_tol + 1e-6)) - 0.5 * (m_consistency / (m_tol + 1e-6))

    slope_obs = obs.get("slope_med", 0.0)
    slope_tmp = tmpl.get("slope_med", 0.0)
    slope_term = -abs(slope_obs - slope_tmp) / (slope_tol + 1e-6)

    knee_obs = obs.get("knee_idx", 0.0)
    knee_tmp = tmpl.get("knee_idx", 0.0)
    knee_term = -abs(knee_obs - knee_tmp) / (knee_tol + 1e-6)

    energy_term = 0.25 * (obs.get("m_ratio", 0.0) + tmpl.get("m_ratio", 0.0))
    ridge_term = 0.1 * (obs.get("ridge_energy", 0.0) + tmpl.get("ridge_energy", 0.0))
    return float(mode_term + slope_term + knee_term + energy_term + ridge_term)


def baseline_clahe(img: np.ndarray) -> np.ndarray:
    clahe = cv2.createCLAHE(clipLimit=2.0, tileGridSize=(8, 8))
    enhanced = clahe.apply((img * 255).astype(np.uint8))
    return normalize_image(enhanced.astype(np.float32))


def baseline_unsharp(img: np.ndarray, amount: float = 1.0, sigma: float = 1.0) -> np.ndarray:
    blurred = cv2.GaussianBlur(img, (0, 0), sigma)
    sharp = cv2.addWeighted(img, 1 + amount, blurred, -amount, 0)
    return normalize_image(sharp)


# --------------------------
# Experiment runner
# --------------------------
def candidate_grid(ells: Iterable[int], pitches: Iterable[float], knees: Iterable[float]) -> List[MaskParams]:
    grid: List[MaskParams] = []
    for ell in ells:
        if ell == 0:
            continue  # avoid trivial ell=0 collapse
        for pitch in pitches:
            for knee in knees:
                grid.append(MaskParams(kind="spiral", ell=ell, chirality=1, pitch=pitch, knee=knee))
                grid.append(MaskParams(kind="spiral", ell=ell, chirality=-1, pitch=pitch, knee=knee))
    # Keep random masks only as baselines; they are not used for selection scoring.
    return grid


def evaluate_decodes(
    ref: np.ndarray, decoded: Dict[str, np.ndarray]
) -> Dict[str, float]:
    metrics: Dict[str, float] = {}
    for name, img in decoded.items():
        metrics[f"{name}_ssim"] = ssim_simple(ref, img)
        metrics[f"{name}_psnr"] = psnr(ref, img)
        metrics[f"{name}_ncc"] = ncc(ref, img)
    return metrics


def run_trial(
    pattern_name: str,
    pattern_img: np.ndarray,
    true_mask: MaskParams,
    decoy_mask: MaskParams,
    noise: NoiseParams,
    candidates: Sequence[MaskParams],
    rng: np.random.Generator,
    ref_amp: float,
    ref_kx: float,
    ref_ky: float,
    approach: str,
    ang_weight: float,
    use_m_distance: bool,
    physics_weight: float,
) -> Dict[str, object]:
    """
    Encode with a true spiral mask, mix in a decoy encoding, decode blindly, and score.
    """
    true_scrambled_1 = encode_image(pattern_img, true_mask, noise, rng, ref_amp=ref_amp, ref_kx=ref_kx, ref_ky=ref_ky)
    true_scrambled_2 = encode_image(pattern_img, true_mask, noise, rng, ref_amp=ref_amp, ref_kx=ref_kx, ref_ky=ref_ky)
    decoy_scrambled_1 = encode_image(pattern_img, decoy_mask, noise, rng, ref_amp=ref_amp, ref_kx=ref_kx, ref_ky=ref_ky)
    decoy_scrambled_2 = encode_image(pattern_img, decoy_mask, noise, rng, ref_amp=ref_amp, ref_kx=ref_kx, ref_ky=ref_ky)

    center = estimate_center(true_scrambled_1)
    lp_n_r, lp_n_theta = 256, 512
    obs_bands_true = annular_band_peaks(log_polar_warp(true_scrambled_1, center=center, n_r=lp_n_r, n_theta=lp_n_theta), n_bands=6)
    obs_bands_decoy = annular_band_peaks(log_polar_warp(decoy_scrambled_1, center=center, n_r=lp_n_r, n_theta=lp_n_theta), n_bands=6)
    obs_features_true = spiral_logpolar_features(true_scrambled_1, center=center, n_r=lp_n_r, n_theta=lp_n_theta)
    obs_features_decoy = spiral_logpolar_features(decoy_scrambled_1, center=center, n_r=lp_n_r, n_theta=lp_n_theta)
    template_cache: Dict[Tuple[float, ...], Dict[str, float]] = {}
    probe_img = template_probe(pattern_img.shape)
    sideband_diag_true: Dict[str, float] = {}
    sideband_diag_decoy: Dict[str, float] = {}
    winding_true = float("nan")
    winding_decoy = float("nan")
    circ_true: Dict[str, float] = {}
    circ_decoy: Dict[str, float] = {}
    charge_true: Dict[str, float] = {}
    charge_decoy: Dict[str, float] = {}
    oam_pre_true = {}
    oam_pre_decoy = {}
    oam_post_true = {}
    oam_post_decoy = {}
    oam_proj_true: Optional[Dict[str, float]] = None
    oam_proj_decoy: Optional[Dict[str, float]] = None
    proj_true: Dict[str, float] = {}
    proj_decoy: Dict[str, float] = {}
    spec_true: Dict[str, float] = {}
    spec_decoy: Dict[str, float] = {}
    knee_amp_true = float("nan")
    knee_amp_decoy = float("nan")
    knee_phasegrad_true = float("nan")
    knee_phasegrad_decoy = float("nan")
    leak_penalty = 0.0
    leak_max = 0.2
    leak_fail = False
    if ref_amp > 0.0:
        def compute_oam_with_gate(field: np.ndarray, center: Tuple[float, float]) -> Dict[str, float]:
            raw = oam_spectrum_annulus(
                field,
                center=center,
                r_inner_frac=0.08,
                r_outer_frac=0.4,
                n_rings=16,
                n_theta=256,
                m_max=3,
            )
            raw_peak = raw.get("m_peak", float("nan"))
            raw_ratio = raw.get("m_peak_ratio", float("nan"))
            total_power = raw.get("total_power", 0.0)
            gate_reason = ""
            gated = False
            m_peak = raw_peak
            m_ratio = raw_ratio
            if not np.isfinite(raw_peak) or not np.isfinite(raw_ratio):
                gated = True
                gate_reason = "nan"
                m_peak = raw_peak if np.isfinite(raw_peak) else 0.0
                m_ratio = raw_ratio if np.isfinite(raw_ratio) else 0.0
            elif total_power <= 1e-8:
                gated = True
                gate_reason = "low_power"
                m_peak = 0.0
                m_ratio = 0.0
            elif raw_ratio < 0.01:
                gated = True
                gate_reason = "flat"
                m_peak = 0.0
                m_ratio = 0.0
            if gate_reason == "":
                gate_reason = "none"
            return {
                "m_peak": m_peak,
                "m_peak_ratio": m_ratio,
                "m_peak_raw": raw_peak,
                "m_peak_ratio_raw": raw_ratio,
                "m_peak_power": raw.get("m_peak_power", float("nan")),
                "total_power": total_power,
                "gated": gated,
                "gate_reason": gate_reason,
            }

        def compute_circulation(field: np.ndarray, center: Tuple[float, float]) -> Dict[str, float]:
            circ = phase_circulation_annulus(
                field,
                center=center,
                r_inner_frac=0.08,
                r_outer_frac=0.4,
                n_rings=16,
                n_theta=256,
            )
            return circ

        def decode_fn(img: np.ndarray, mp: MaskParams):
            dec, _ = decode_offaxis(img, mp, ref_kx=ref_kx, ref_ky=ref_ky, auto_sideband=True)
            return dec
        field_true, diag_true = reconstruct_complex_field(true_scrambled_1, ref_kx=ref_kx, ref_ky=ref_ky, auto_sideband=True)
        field_true = remove_phase_ramp(field_true)
        phase_std_true, imag_max_true = ensure_phase_alive(field_true, tag="charge_tap_true")
        fp_field_true = field_fingerprint(field_true)
        core_true = estimate_vortex_center(field_true)
        charge_true = charge_metrics(field_true, inner_margin_frac=0.05, center=core_true)
        charge_true["charge_phase_std"] = charge_true.get("phase_std", float("nan"))
        charge_true["charge_imag_max"] = charge_true.get("imag_max", float("nan"))
        charge_true["charge_phase_std_raw"] = phase_std_true
        charge_true["charge_imag_max_raw"] = imag_max_true
        diag_true.update(
            {
                "vortex_cx": core_true[0],
                "vortex_cy": core_true[1],
                "charge_sum": charge_true.get("sum_charge", float("nan")),
                "charge_sum_abs": charge_true.get("sum_charge_abs", float("nan")),
                "charge_n_vort": charge_true.get("n_vort", float("nan")),
                "charge_cx": charge_true.get("charge_cx", float("nan")),
                "charge_cy": charge_true.get("charge_cy", float("nan")),
                "q_float_mean": charge_true.get("q_float_mean", float("nan")),
                "q_float_std": charge_true.get("q_float_std", float("nan")),
                "q_float_max_abs": charge_true.get("q_float_max_abs", float("nan")),
                "n_q_gt_0p25": charge_true.get("n_q_gt_0p25", float("nan")),
                "n_q_gt_0p5": charge_true.get("n_q_gt_0p5", float("nan")),
                "q_abs_p50": charge_true.get("q_abs_p50", float("nan")),
                "q_abs_p90": charge_true.get("q_abs_p90", float("nan")),
                "n_valid_cells": charge_true.get("n_valid_cells", float("nan")),
                "valid_cell_frac": charge_true.get("valid_cell_frac", float("nan")),
                "charge_phase_std": charge_true.get("phase_std", float("nan")),
                "charge_imag_max": charge_true.get("imag_max", float("nan")),
                # raw tripwire stats before charge mapping (same field)
                "charge_phase_std_raw": phase_std_true,
                "charge_imag_max_raw": imag_max_true,
                # fingerprints to catch field swaps
                "fp_mean_real": fp_field_true["mean_real"],
                "fp_std_real": fp_field_true["std_real"],
                "fp_mean_imag": fp_field_true["mean_imag"],
                "fp_std_imag": fp_field_true["std_imag"],
                "fp_mean_abs": fp_field_true["mean_abs"],
            }
        )
        # OAM around charge centroid to be robust to splitting
        # Fallback priority for OAM center: significant charge centroid -> vortex core -> image center
        oam_cx_true = charge_true.get("charge_cx", float("nan"))
        oam_cy_true = charge_true.get("charge_cy", float("nan"))
        if not np.isfinite(oam_cx_true) or not np.isfinite(oam_cy_true):
            oam_cx_true, oam_cy_true = core_true
        oam_center_true = (oam_cx_true, oam_cy_true)
        proj_scan_offsets = (
            (0.0, 0.0),
            (1.0, 0.0),
            (-1.0, 0.0),
            (0.0, 1.0),
            (0.0, -1.0),
            (1.0, 1.0),
            (1.0, -1.0),
            (-1.0, 1.0),
            (-1.0, -1.0),
        )
        proj_true = helical_projection(
            field_true,
            center=oam_center_true,
            ells=(-3, -2, -1, 1, 2, 3),
            r_inner_frac=0.08,
            r_outer_frac=0.4,
            p=0.25,
            n_rings=3,
            phase_only=True,
            scan_offsets=proj_scan_offsets,
        )
        diag_true.update(proj_true)
        spec_true = helical_spectrum_stability(
            field_true,
            center=oam_center_true,
            ells=(-3, -2, -1, 1, 2, 3),
        )
        diag_true.update(spec_true)
        sideband_diag_true = diag_true
        winding_true = phase_winding(field_true, center=core_true)
        oam_pre_true = compute_oam_with_gate(field_true, center=oam_center_true)
        circ_true = compute_circulation(field_true, center=oam_center_true)
        # Promote the projection witness to an OAM-like target: pick the ell with the
        # largest projection magnitude and use a peak/second ratio as reliability.
        def proj_to_oam_target(proj: Dict[str, float]) -> Dict[str, float]:
            keys = [k for k in proj.keys() if k.startswith("proj_ell_")]
            vals = []
            for k in keys:
                try:
                    ell = int(k.split("proj_ell_")[-1])
                except Exception:
                    continue
                vals.append((ell, float(proj[k])))
            if not vals:
                return {
                    "m_peak": 0.0,
                    "m_peak_ratio": 0.0,
                    "proj_peak_ratio": 0.0,
                    "gated": True,
                    "gate_reason": "no_proj",
                }
            vals.sort(key=lambda t: t[1], reverse=True)
            peak_ell, peak_val = vals[0]
            second_val = vals[1][1] if len(vals) > 1 else 1e-12
            ratio = peak_val / (second_val + 1e-9)
            return {
                "m_peak": float(peak_ell),
                "m_peak_ratio": float(ratio),
                "m_peak_raw": float(peak_ell),
                "m_peak_ratio_raw": float(ratio),
                "proj_peak_ell": float(peak_ell),
                "proj_peak_val": float(peak_val),
                "proj_second_val": float(second_val),
                "proj_peak_ratio": float(ratio),
                "proj_peak_margin": float(peak_val - second_val),
                "gated": False,
                "gate_reason": "proj",
            }

        oam_proj_true = proj_to_oam_target(proj_true)
        field_decoy, diag_decoy = reconstruct_complex_field(decoy_scrambled_1, ref_kx=ref_kx, ref_ky=ref_ky, auto_sideband=True)
        field_decoy = remove_phase_ramp(field_decoy)
        phase_std_decoy, imag_max_decoy = ensure_phase_alive(field_decoy, tag="charge_tap_decoy")
        fp_field_decoy = field_fingerprint(field_decoy)
        core_decoy = estimate_vortex_center(field_decoy)
        charge_decoy = charge_metrics(field_decoy, inner_margin_frac=0.05, center=core_decoy)
        charge_decoy["charge_phase_std"] = charge_decoy.get("phase_std", float("nan"))
        charge_decoy["charge_imag_max"] = charge_decoy.get("imag_max", float("nan"))
        charge_decoy["charge_phase_std_raw"] = phase_std_decoy
        charge_decoy["charge_imag_max_raw"] = imag_max_decoy
        diag_decoy.update(
            {
                "vortex_cx": core_decoy[0],
                "vortex_cy": core_decoy[1],
                "charge_sum": charge_decoy.get("sum_charge", float("nan")),
                "charge_sum_abs": charge_decoy.get("sum_charge_abs", float("nan")),
                "charge_n_vort": charge_decoy.get("n_vort", float("nan")),
                "charge_cx": charge_decoy.get("charge_cx", float("nan")),
                "charge_cy": charge_decoy.get("charge_cy", float("nan")),
                "q_float_mean": charge_decoy.get("q_float_mean", float("nan")),
                "q_float_std": charge_decoy.get("q_float_std", float("nan")),
                "q_float_max_abs": charge_decoy.get("q_float_max_abs", float("nan")),
                "n_q_gt_0p25": charge_decoy.get("n_q_gt_0p25", float("nan")),
                "n_q_gt_0p5": charge_decoy.get("n_q_gt_0p5", float("nan")),
                "q_abs_p50": charge_decoy.get("q_abs_p50", float("nan")),
                "q_abs_p90": charge_decoy.get("q_abs_p90", float("nan")),
                "n_valid_cells": charge_decoy.get("n_valid_cells", float("nan")),
                "valid_cell_frac": charge_decoy.get("valid_cell_frac", float("nan")),
                "charge_phase_std": charge_decoy.get("phase_std", float("nan")),
                "charge_imag_max": charge_decoy.get("imag_max", float("nan")),
                "charge_phase_std_raw": phase_std_decoy,
                "charge_imag_max_raw": imag_max_decoy,
                "fp_mean_real": fp_field_decoy["mean_real"],
                "fp_std_real": fp_field_decoy["std_real"],
                "fp_mean_imag": fp_field_decoy["mean_imag"],
                "fp_std_imag": fp_field_decoy["std_imag"],
                "fp_mean_abs": fp_field_decoy["mean_abs"],
            }
        )
        # OAM center for decoy: prefer charge centroid, then vortex core
        oam_cx_decoy = charge_decoy.get("charge_cx", float("nan"))
        oam_cy_decoy = charge_decoy.get("charge_cy", float("nan"))
        if not np.isfinite(oam_cx_decoy) or not np.isfinite(oam_cy_decoy):
            oam_cx_decoy, oam_cy_decoy = core_decoy
        oam_center_decoy = (oam_cx_decoy, oam_cy_decoy)
        proj_decoy = helical_projection(
            field_decoy,
            center=oam_center_decoy,
            ells=(-3, -2, -1, 1, 2, 3),
            r_inner_frac=0.08,
            r_outer_frac=0.4,
            p=0.25,
            n_rings=3,
            phase_only=True,
            scan_offsets=proj_scan_offsets,
        )
        diag_decoy.update(proj_decoy)
        spec_decoy = helical_spectrum_stability(
            field_decoy,
            center=oam_center_decoy,
            ells=(-3, -2, -1, 1, 2, 3),
        )
        diag_decoy.update(spec_decoy)
        sideband_diag_decoy = diag_decoy
        winding_decoy = phase_winding(field_decoy, center=core_decoy)
        oam_pre_decoy = compute_oam_with_gate(field_decoy, center=oam_center_decoy)
        circ_decoy = compute_circulation(field_decoy, center=oam_center_decoy)
        oam_proj_decoy = proj_to_oam_target(proj_decoy)
        leak_val = max(diag_true.get("twin_leak", 0.0), diag_decoy.get("twin_leak", 0.0))
        leak_penalty = max(0.0, (leak_val - leak_max) / max(leak_max, 1e-6))
        leak_fail = leak_penalty > 0.0
        amp_true = np.abs(field_true).astype(np.float32)
        amp_decoy = np.abs(field_decoy).astype(np.float32)
        amp_lp_true = log_polar_warp(amp_true, center=center, n_r=lp_n_r, n_theta=lp_n_theta)
        amp_lp_decoy = log_polar_warp(amp_decoy, center=center, n_r=lp_n_r, n_theta=lp_n_theta)
        knee_amp_true = knee_features(amp_lp_true).get("knee_idx", float("nan"))
        knee_amp_decoy = knee_features(amp_lp_decoy).get("knee_idx", float("nan"))
        phase_true = np.angle(field_true).astype(np.float32)
        phase_decoy = np.angle(field_decoy).astype(np.float32)
        gx_t = cv2.Sobel(phase_true, cv2.CV_32F, 1, 0, ksize=3)
        gy_t = cv2.Sobel(phase_true, cv2.CV_32F, 0, 1, ksize=3)
        phase_grad_true = np.sqrt(gx_t ** 2 + gy_t ** 2)
        gx_d = cv2.Sobel(phase_decoy, cv2.CV_32F, 1, 0, ksize=3)
        gy_d = cv2.Sobel(phase_decoy, cv2.CV_32F, 0, 1, ksize=3)
        phase_grad_decoy = np.sqrt(gx_d ** 2 + gy_d ** 2)
        pg_lp_true = log_polar_warp(phase_grad_true, center=center, n_r=lp_n_r, n_theta=lp_n_theta)
        pg_lp_decoy = log_polar_warp(phase_grad_decoy, center=center, n_r=lp_n_r, n_theta=lp_n_theta)
        knee_phasegrad_true = knee_features(pg_lp_true).get("knee_idx", float("nan"))
        knee_phasegrad_decoy = knee_features(pg_lp_decoy).get("knee_idx", float("nan"))
    else:
        decode_fn = decode_with_mask

    decoded_true, best_mask_true, self_score_true, margin_true, top_true, physics_score_true = grid_decode_joint(
        [true_scrambled_1, true_scrambled_2],
        candidates,
        center=center,
        n_r=lp_n_r,
        n_theta=lp_n_theta,
        obs_bands=obs_bands_true,
        ref_amp=ref_amp,
        ref_kx=ref_kx,
        ref_ky=ref_ky,
        ang_weight=ang_weight,
        use_m_distance=use_m_distance,
        obs_features=obs_features_true,
        physics_weight=physics_weight,
        template_cache=template_cache,
        probe=probe_img,
        decode_fn=decode_fn,
        leak_penalty=leak_penalty,
        target_winding=oam_proj_true.get("m_peak") if isinstance(oam_proj_true, dict) else None,
        # merge projection witness into oam_target so the scoring term can see proj_ell_* fields
        oam_target={**oam_proj_true, **proj_true} if isinstance(oam_proj_true, dict) else None,
        oam_weight=1.0,
    )
    decoded_decoy, best_mask_decoy, self_score_decoy, margin_decoy, top_decoy, physics_score_decoy = grid_decode_joint(
        [decoy_scrambled_1, decoy_scrambled_2],
        candidates,
        center=center,
        n_r=lp_n_r,
        n_theta=lp_n_theta,
        obs_bands=obs_bands_decoy,
        ref_amp=ref_amp,
        ref_kx=ref_kx,
        ref_ky=ref_ky,
        ang_weight=ang_weight,
        use_m_distance=use_m_distance,
        obs_features=obs_features_decoy,
        physics_weight=physics_weight,
        template_cache=template_cache,
        probe=probe_img,
        decode_fn=decode_fn,
        leak_penalty=leak_penalty,
        target_winding=oam_proj_decoy.get("m_peak") if isinstance(oam_proj_decoy, dict) else None,
        oam_target={**oam_proj_decoy, **proj_decoy} if isinstance(oam_proj_decoy, dict) else None,
        oam_weight=1.0,
    )

    baselines_true = {
        "clahe": baseline_clahe(true_scrambled_1),
        "unsharp": baseline_unsharp(true_scrambled_1),
    }
    baselines_decoy = {
        "clahe": baseline_clahe(decoy_scrambled_1),
        "unsharp": baseline_unsharp(decoy_scrambled_1),
    }

    # Log-polar annular coherence (use first scramble for each)
    center = (pattern_img.shape[1] / 2.0, pattern_img.shape[0] / 2.0)
    lp_true = log_polar_warp(true_scrambled_1, center=center, n_r=256, n_theta=256)
    lp_decoy = log_polar_warp(decoy_scrambled_1, center=center, n_r=256, n_theta=256)
    bands_true = annular_band_peaks(lp_true, n_bands=6)
    bands_decoy = annular_band_peaks(lp_decoy, n_bands=6)
    band_summary_true = summarize_bands(bands_true)
    band_summary_decoy = summarize_bands(bands_decoy)

    metrics_true = evaluate_decodes(pattern_img, {"matched": decoded_true, **baselines_true})
    metrics_decoy = evaluate_decodes(pattern_img, {"matched": decoded_decoy, **baselines_decoy})
    # Post-decode OAM (diagnostic only; should collapse toward 0 for correct masks)
    try:
        oam_post_true = oam_spectrum(decoded_true.astype(np.complex64), center=center, n_r=128, n_theta=256, m_max=8)
    except Exception:
        oam_post_true = {}
    try:
        oam_post_decoy = oam_spectrum(decoded_decoy.astype(np.complex64), center=center, n_r=128, n_theta=256, m_max=8)
    except Exception:
        oam_post_decoy = {}

    # Diagnostic: likelihood rank of true mask using known ground truth (ablation)
    true_like_scores = []
    for params in candidates:
        reenc = forward_intensity(pattern_img, params, ref_amp=ref_amp, ref_kx=ref_kx, ref_ky=ref_ky)
        mse = float(np.mean((true_scrambled_1 - reenc) ** 2))
        true_like_scores.append((-mse, params))
    true_like_scores.sort(key=lambda x: x[0], reverse=True)
    true_like_rank = next(
        (idx for idx, (_, p) in enumerate(true_like_scores) if p == true_mask), len(true_like_scores)
    )
    best_like_score = true_like_scores[0][0]
    true_like_score = next((s for s, p in true_like_scores if p == true_mask), float("-inf"))
    like_gap = best_like_score - true_like_score

    return {
        "pattern": pattern_name,
        "true_mask": asdict(true_mask),
        "decoy_mask": asdict(decoy_mask),
        "best_mask_true": asdict(best_mask_true),
        "best_mask_decoy": asdict(best_mask_decoy),
        "self_score_true": self_score_true,
        "self_score_decoy": self_score_decoy,
        "margin_true": margin_true,
        "margin_decoy": margin_decoy,
        "top_scores_true": top_true,
        "top_scores_decoy": top_decoy,
        "true_like_rank": true_like_rank,
        "true_like_gap": like_gap,
        "bands_true": bands_true,
        "bands_decoy": bands_decoy,
        "band_summary_true": band_summary_true,
        "band_summary_decoy": band_summary_decoy,
        "metrics_true": metrics_true,
        "metrics_decoy": metrics_decoy,
        "physics_score_true": physics_score_true,
        "physics_score_decoy": physics_score_decoy,
        "obs_features_true": obs_features_true,
        "obs_features_decoy": obs_features_decoy,
        "proj_true": proj_true if ref_amp > 0.0 else {},
        "proj_decoy": proj_decoy if ref_amp > 0.0 else {},
        "spec_true": spec_true if ref_amp > 0.0 else {},
        "spec_decoy": spec_decoy if ref_amp > 0.0 else {},
        "sideband_diag_true": sideband_diag_true,
        "sideband_diag_decoy": sideband_diag_decoy,
        "winding_true": winding_true,
        "winding_decoy": winding_decoy,
        "circ_true": circ_true if ref_amp > 0.0 else {},
        "circ_decoy": circ_decoy if ref_amp > 0.0 else {},
        "vortex_cx_true": sideband_diag_true.get("vortex_cx", float("nan")) if isinstance(sideband_diag_true, dict) else float("nan"),
        "vortex_cy_true": sideband_diag_true.get("vortex_cy", float("nan")) if isinstance(sideband_diag_true, dict) else float("nan"),
        "vortex_cx_decoy": sideband_diag_decoy.get("vortex_cx", float("nan")) if isinstance(sideband_diag_decoy, dict) else float("nan"),
        "vortex_cy_decoy": sideband_diag_decoy.get("vortex_cy", float("nan")) if isinstance(sideband_diag_decoy, dict) else float("nan"),
        "charge_diag_true": charge_true if ref_amp > 0.0 else {},
        "charge_diag_decoy": charge_decoy if ref_amp > 0.0 else {},
        "knee_amp_true": knee_amp_true,
        "knee_amp_decoy": knee_amp_decoy,
        "knee_phasegrad_true": knee_phasegrad_true,
        "knee_phasegrad_decoy": knee_phasegrad_decoy,
        "leak_penalty": leak_penalty,
        "leak_fail": leak_fail,
        "oam_pre_true": oam_pre_true,
        "oam_pre_decoy": oam_pre_decoy,
        "peak_y_true": sideband_diag_true.get("peak_y", float("nan")) if isinstance(sideband_diag_true, dict) else float("nan"),
        "peak_x_true": sideband_diag_true.get("peak_x", float("nan")) if isinstance(sideband_diag_true, dict) else float("nan"),
        "crop_radius_true": sideband_diag_true.get("crop_radius", float("nan")) if isinstance(sideband_diag_true, dict) else float("nan"),
        "peak_y_decoy": sideband_diag_decoy.get("peak_y", float("nan")) if isinstance(sideband_diag_decoy, dict) else float("nan"),
        "peak_x_decoy": sideband_diag_decoy.get("peak_x", float("nan")) if isinstance(sideband_diag_decoy, dict) else float("nan"),
        "crop_radius_decoy": sideband_diag_decoy.get("crop_radius", float("nan")) if isinstance(sideband_diag_decoy, dict) else float("nan"),
        "peak_y_pred_true": sideband_diag_true.get("peak_y_pred", float("nan")) if isinstance(sideband_diag_true, dict) else float("nan"),
        "peak_x_pred_true": sideband_diag_true.get("peak_x_pred", float("nan")) if isinstance(sideband_diag_true, dict) else float("nan"),
        "peak_y_pred_decoy": sideband_diag_decoy.get("peak_y_pred", float("nan")) if isinstance(sideband_diag_decoy, dict) else float("nan"),
        "peak_x_pred_decoy": sideband_diag_decoy.get("peak_x_pred", float("nan")) if isinstance(sideband_diag_decoy, dict) else float("nan"),
        "r_min_true": sideband_diag_true.get("r_min", float("nan")) if isinstance(sideband_diag_true, dict) else float("nan"),
        "r_min_decoy": sideband_diag_decoy.get("r_min", float("nan")) if isinstance(sideband_diag_decoy, dict) else float("nan"),
        "dist_detect_pred_plus_true": sideband_diag_true.get("dist_detect_pred_plus", float("nan")) if isinstance(sideband_diag_true, dict) else float("nan"),
        "dist_detect_pred_minus_true": sideband_diag_true.get("dist_detect_pred_minus", float("nan")) if isinstance(sideband_diag_true, dict) else float("nan"),
        "dist_detect_pred_plus_decoy": sideband_diag_decoy.get("dist_detect_pred_plus", float("nan")) if isinstance(sideband_diag_decoy, dict) else float("nan"),
        "dist_detect_pred_minus_decoy": sideband_diag_decoy.get("dist_detect_pred_minus", float("nan")) if isinstance(sideband_diag_decoy, dict) else float("nan"),
        "kx_used_bins_true": sideband_diag_true.get("kx_used_bins", float("nan")) if isinstance(sideband_diag_true, dict) else float("nan"),
        "ky_used_bins_true": sideband_diag_true.get("ky_used_bins", float("nan")) if isinstance(sideband_diag_true, dict) else float("nan"),
        "kx_used_bins_decoy": sideband_diag_decoy.get("kx_used_bins", float("nan")) if isinstance(sideband_diag_decoy, dict) else float("nan"),
        "ky_used_bins_decoy": sideband_diag_decoy.get("ky_used_bins", float("nan")) if isinstance(sideband_diag_decoy, dict) else float("nan"),
        "charge_sum_true": charge_true.get("sum_charge", float("nan")),
        "charge_sum_abs_true": charge_true.get("sum_charge_abs", float("nan")),
        "charge_n_vort_true": charge_true.get("n_vort", float("nan")),
        "charge_cx_true": charge_true.get("charge_cx", float("nan")),
        "charge_cy_true": charge_true.get("charge_cy", float("nan")),
        "q_float_mean_true": charge_true.get("q_float_mean", float("nan")),
        "q_float_std_true": charge_true.get("q_float_std", float("nan")),
        "q_float_max_abs_true": charge_true.get("q_float_max_abs", float("nan")),
        "n_q_gt_0p25_true": charge_true.get("n_q_gt_0p25", float("nan")),
        "n_q_gt_0p5_true": charge_true.get("n_q_gt_0p5", float("nan")),
        "q_abs_p50_true": charge_true.get("q_abs_p50", float("nan")),
        "q_abs_p90_true": charge_true.get("q_abs_p90", float("nan")),
        "n_valid_cells_true": charge_true.get("n_valid_cells", float("nan")),
        "valid_cell_frac_true": charge_true.get("valid_cell_frac", float("nan")),
        "fp_mean_real_true": sideband_diag_true.get("fp_mean_real", float("nan")),
        "fp_std_real_true": sideband_diag_true.get("fp_std_real", float("nan")),
        "fp_mean_imag_true": sideband_diag_true.get("fp_mean_imag", float("nan")),
        "fp_std_imag_true": sideband_diag_true.get("fp_std_imag", float("nan")),
        "fp_mean_abs_true": sideband_diag_true.get("fp_mean_abs", float("nan")),
        "charge_sum_decoy": charge_decoy.get("sum_charge", float("nan")),
        "charge_sum_abs_decoy": charge_decoy.get("sum_charge_abs", float("nan")),
        "charge_n_vort_decoy": charge_decoy.get("n_vort", float("nan")),
        "charge_cx_decoy": charge_decoy.get("charge_cx", float("nan")),
        "charge_cy_decoy": charge_decoy.get("charge_cy", float("nan")),
        "q_float_mean_decoy": charge_decoy.get("q_float_mean", float("nan")),
        "q_float_std_decoy": charge_decoy.get("q_float_std", float("nan")),
        "q_float_max_abs_decoy": charge_decoy.get("q_float_max_abs", float("nan")),
        "n_q_gt_0p25_decoy": charge_decoy.get("n_q_gt_0p25", float("nan")),
        "n_q_gt_0p5_decoy": charge_decoy.get("n_q_gt_0p5", float("nan")),
        "q_abs_p50_decoy": charge_decoy.get("q_abs_p50", float("nan")),
        "q_abs_p90_decoy": charge_decoy.get("q_abs_p90", float("nan")),
        "n_valid_cells_decoy": charge_decoy.get("n_valid_cells", float("nan")),
        "valid_cell_frac_decoy": charge_decoy.get("valid_cell_frac", float("nan")),
        "fp_mean_real_decoy": sideband_diag_decoy.get("fp_mean_real", float("nan")),
        "fp_std_real_decoy": sideband_diag_decoy.get("fp_std_real", float("nan")),
        "fp_mean_imag_decoy": sideband_diag_decoy.get("fp_mean_imag", float("nan")),
        "fp_std_imag_decoy": sideband_diag_decoy.get("fp_std_imag", float("nan")),
        "fp_mean_abs_decoy": sideband_diag_decoy.get("fp_mean_abs", float("nan")),
        "oam_peak_pre_true": oam_pre_true.get("m_peak", float("nan")),
        "oam_peak_ratio_pre_true": oam_pre_true.get("m_peak_ratio", float("nan")),
        "oam_pre_raw_peak_true": oam_pre_true.get("m_peak_raw", float("nan")),
        "oam_pre_raw_ratio_true": oam_pre_true.get("m_peak_ratio_raw", float("nan")),
        "oam_pre_gated_true": int(bool(oam_pre_true.get("gated", False))),
        "oam_pre_gate_reason_true": oam_pre_true.get("gate_reason", ""),
        "oam_peak_pre_decoy": oam_pre_decoy.get("m_peak", float("nan")),
        "oam_peak_ratio_pre_decoy": oam_pre_decoy.get("m_peak_ratio", float("nan")),
        "oam_pre_raw_peak_decoy": oam_pre_decoy.get("m_peak_raw", float("nan")),
        "oam_pre_raw_ratio_decoy": oam_pre_decoy.get("m_peak_ratio_raw", float("nan")),
        "oam_pre_gated_decoy": int(bool(oam_pre_decoy.get("gated", False))),
        "oam_pre_gate_reason_decoy": oam_pre_decoy.get("gate_reason", ""),
        "oam_peak_post_true": oam_post_true.get("m_peak", float("nan")) if isinstance(oam_post_true, dict) else float("nan"),
        "oam_peak_ratio_post_true": oam_post_true.get("m_peak_ratio", float("nan")) if isinstance(oam_post_true, dict) else float("nan"),
        "oam_peak_post_decoy": oam_post_decoy.get("m_peak", float("nan")) if isinstance(oam_post_decoy, dict) else float("nan"),
        "oam_peak_ratio_post_decoy": oam_post_decoy.get("m_peak_ratio", float("nan")) if isinstance(oam_post_decoy, dict) else float("nan"),
        "circ_m_med_true": circ_true.get("m_circ_med", float("nan")) if circ_true else float("nan"),
        "circ_m_std_true": circ_true.get("m_circ_std", float("nan")) if circ_true else float("nan"),
        "circ_m_n_true": circ_true.get("m_circ_n", float("nan")) if circ_true else float("nan"),
        "m_circ_kept_frac_true": circ_true.get("m_circ_kept_frac", float("nan")) if circ_true else float("nan"),
        "circ_m_med_decoy": circ_decoy.get("m_circ_med", float("nan")) if circ_decoy else float("nan"),
        "circ_m_std_decoy": circ_decoy.get("m_circ_std", float("nan")) if circ_decoy else float("nan"),
        "circ_m_n_decoy": circ_decoy.get("m_circ_n", float("nan")) if circ_decoy else float("nan"),
        "m_circ_kept_frac_decoy": circ_decoy.get("m_circ_kept_frac", float("nan")) if circ_decoy else float("nan"),
        # Flatten projection witness so CSV columns exist
        "proj_ell_-3_true": proj_true.get("proj_ell_-3", float("nan")),
        "proj_ell_-2_true": proj_true.get("proj_ell_-2", float("nan")),
        "proj_ell_-1_true": proj_true.get("proj_ell_-1", float("nan")),
        "proj_ell_1_true": proj_true.get("proj_ell_1", float("nan")),
        "proj_ell_2_true": proj_true.get("proj_ell_2", float("nan")),
        "proj_ell_3_true": proj_true.get("proj_ell_3", float("nan")),
        "proj_ell_-3_decoy": proj_decoy.get("proj_ell_-3", float("nan")),
        "proj_ell_-2_decoy": proj_decoy.get("proj_ell_-2", float("nan")),
        "proj_ell_-1_decoy": proj_decoy.get("proj_ell_-1", float("nan")),
        "proj_ell_1_decoy": proj_decoy.get("proj_ell_1", float("nan")),
        "proj_ell_2_decoy": proj_decoy.get("proj_ell_2", float("nan")),
        "proj_ell_3_decoy": proj_decoy.get("proj_ell_3", float("nan")),
        "proj_peak_true": oam_proj_true.get("proj_peak_ell", float("nan")) if isinstance(oam_proj_true, dict) else float("nan"),
        "proj_peak_ratio_true": oam_proj_true.get("proj_peak_ratio", float("nan")) if isinstance(oam_proj_true, dict) else float("nan"),
        "proj_peak_decoy": oam_proj_decoy.get("proj_peak_ell", float("nan")) if isinstance(oam_proj_decoy, dict) else float("nan"),
        "proj_peak_ratio_decoy": oam_proj_decoy.get("proj_peak_ratio", float("nan")) if isinstance(oam_proj_decoy, dict) else float("nan"),
        "proj_center_dx_true": proj_true.get("proj_center_dx", float("nan")),
        "proj_center_dy_true": proj_true.get("proj_center_dy", float("nan")),
        "proj_center_dx_decoy": proj_decoy.get("proj_center_dx", float("nan")),
        "proj_center_dy_decoy": proj_decoy.get("proj_center_dy", float("nan")),
        "spec_peak_ell_true": spec_true.get("spec_peak_ell", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_peak_ratio_true": spec_true.get("spec_peak_ratio", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_stability_true": spec_true.get("spec_stability_frac", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_peak_ell_decoy": spec_decoy.get("spec_peak_ell", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_peak_ratio_decoy": spec_decoy.get("spec_peak_ratio", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_stability_decoy": spec_decoy.get("spec_stability_frac", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_med_ell_-3_true": spec_true.get("spec_med_ell_-3", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_med_ell_-2_true": spec_true.get("spec_med_ell_-2", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_med_ell_-1_true": spec_true.get("spec_med_ell_-1", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_med_ell_1_true": spec_true.get("spec_med_ell_1", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_med_ell_2_true": spec_true.get("spec_med_ell_2", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_med_ell_3_true": spec_true.get("spec_med_ell_3", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_std_ell_-3_true": spec_true.get("spec_std_ell_-3", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_std_ell_-2_true": spec_true.get("spec_std_ell_-2", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_std_ell_-1_true": spec_true.get("spec_std_ell_-1", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_std_ell_1_true": spec_true.get("spec_std_ell_1", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_std_ell_2_true": spec_true.get("spec_std_ell_2", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_std_ell_3_true": spec_true.get("spec_std_ell_3", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_med_ell_-3_decoy": spec_decoy.get("spec_med_ell_-3", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_med_ell_-2_decoy": spec_decoy.get("spec_med_ell_-2", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_med_ell_-1_decoy": spec_decoy.get("spec_med_ell_-1", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_med_ell_1_decoy": spec_decoy.get("spec_med_ell_1", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_med_ell_2_decoy": spec_decoy.get("spec_med_ell_2", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_med_ell_3_decoy": spec_decoy.get("spec_med_ell_3", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_std_ell_-3_decoy": spec_decoy.get("spec_std_ell_-3", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_std_ell_-2_decoy": spec_decoy.get("spec_std_ell_-2", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_std_ell_-1_decoy": spec_decoy.get("spec_std_ell_-1", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_std_ell_1_decoy": spec_decoy.get("spec_std_ell_1", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_std_ell_2_decoy": spec_decoy.get("spec_std_ell_2", float("nan")) if ref_amp > 0.0 else float("nan"),
        "spec_std_ell_3_decoy": spec_decoy.get("spec_std_ell_3", float("nan")) if ref_amp > 0.0 else float("nan"),
        "true_scrambled": true_scrambled_1,
        "decoy_scrambled": decoy_scrambled_1,
        "decoded_true": decoded_true,
        "decoded_decoy": decoded_decoy,
        "approach": approach,
    }


def run_experiment(
    out_dir: Path,
    trials: int,
    image_size: int,
    noise: NoiseParams,
    save_images: bool,
    seed: int,
    pattern_filter: Optional[set] = None,
    ref_amp: float = 0.0,
    ref_kx: float = 0.0,
    ref_ky: float = 0.0,
    approaches: Sequence[str] = ("baseline",),
    ells: Optional[Sequence[int]] = None,
    pitches: Optional[Sequence[float]] = None,
    knees: Optional[Sequence[float]] = None,
) -> pd.DataFrame:
    out_dir.mkdir(parents=True, exist_ok=True)
    date_str = datetime.now().strftime("%y%m%d")
    existing = sorted(out_dir.glob(f"spiral_decoder_results_{date_str}_*.csv"))
    next_run = 1
    if existing:
        try:
            last = existing[-1].stem.split("_")[-1]
            next_run = int(last) + 1
        except Exception:
            next_run = len(existing) + 1
    stamped = out_dir / f"spiral_decoder_results_{date_str}_{next_run:03d}.csv"
    csv_path = out_dir / "spiral_decoder_results.csv"
    images_dir = out_dir / f"run_{date_str}_{next_run:03d}"
    if save_images:
        images_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    patterns = generate_patterns(size=image_size)
    if pattern_filter:
        patterns = {k: v for k, v in patterns.items() if k in pattern_filter}
    grid_ells = [-3, -2, -1, 1, 2, 3] if ells is None else [e for e in ells if e != 0]
    grid_pitches = [0.0, 1.0, 1.5, 2.5, 3.0] if pitches is None else list(pitches)
    grid_knees = [0.0, 0.5, 1.0] if knees is None else list(knees)
    grid = candidate_grid(ells=grid_ells, pitches=grid_pitches, knees=grid_knees)

    rows: List[Dict[str, object]] = []
    for trial in range(trials):
        for name, img in patterns.items():
            # Draw true mask parameters
            ell = rng.integers(-3, 4)
            if ell == 0:
                ell = 1
            true_mask = MaskParams(
                kind="spiral",
                ell=int(ell),
                chirality=int(rng.choice([-1, 1])),
                pitch=float(rng.choice([0.0, 1.0, 2.5])),
                knee=float(rng.choice([0.0, 0.5, 1.0])),
            )
            decoy_mask = MaskParams(
                kind="spiral",
                ell=int(-true_mask.ell),
                chirality=-true_mask.chirality,
                pitch=true_mask.pitch + 1.0,
                knee=true_mask.knee + 0.5,
            )

            for approach in approaches:
                # Approach presets
                approach_lower = approach.strip().lower()
                ref_amp_use = ref_amp
                ang_weight = 0.1
                use_m_distance = False
                physics_weight = 0.5
                if approach_lower == "interferometric":
                    if ref_amp_use <= 0:
                        ref_amp_use = 0.2
                    ang_weight = 0.1
                    use_m_distance = False
                    physics_weight = 0.6
                elif approach_lower == "angular":
                    ang_weight = 0.2
                    use_m_distance = True
                    physics_weight = 0.7
                else:
                    approach_lower = "baseline"

                result = run_trial(
                    pattern_name=name,
                    pattern_img=img,
                    true_mask=true_mask,
                    decoy_mask=decoy_mask,
                    noise=noise,
                    candidates=grid,
                    rng=rng,
                    ref_amp=ref_amp_use,
                    ref_kx=ref_kx,
                    ref_ky=ref_ky,
                    approach=approach_lower,
                    ang_weight=ang_weight,
                    use_m_distance=use_m_distance,
                    physics_weight=physics_weight,
                )

                # Metrics flattening for CSV
                row = {
                    "approach": approach_lower,
                    "trial": trial,
                    "pattern": name,
                    "true_mask": json.dumps(result["true_mask"]),
                    "decoy_mask": json.dumps(result["decoy_mask"]),
                    "best_mask_true": json.dumps(result["best_mask_true"]),
                    "best_mask_decoy": json.dumps(result["best_mask_decoy"]),
                    "self_score_true": result["self_score_true"],
                    "self_score_decoy": result["self_score_decoy"],
                    "margin_true": result["margin_true"],
                    "margin_decoy": result["margin_decoy"],
                    "top_scores_true": json.dumps(result["top_scores_true"]),
                    "top_scores_decoy": json.dumps(result["top_scores_decoy"]),
                    "true_like_rank": result["true_like_rank"],
                    "true_like_gap": result["true_like_gap"],
                    "bands_true": json.dumps(result["bands_true"]),
                    "bands_decoy": json.dumps(result["bands_decoy"]),
                    "band_summary_true": json.dumps(result["band_summary_true"]),
                    "band_summary_decoy": json.dumps(result["band_summary_decoy"]),
                    "physics_score_true": result["physics_score_true"],
                    "physics_score_decoy": result["physics_score_decoy"],
                    "obs_features_true": json.dumps(result["obs_features_true"]),
                    "obs_features_decoy": json.dumps(result["obs_features_decoy"]),
                    "sideband_snr_true": result["sideband_diag_true"].get("sideband_snr", float("nan")) if isinstance(result["sideband_diag_true"], dict) else float("nan"),
                    "sideband_snr_decoy": result["sideband_diag_decoy"].get("sideband_snr", float("nan")) if isinstance(result["sideband_diag_decoy"], dict) else float("nan"),
                    "twin_leak_true": result["sideband_diag_true"].get("twin_leak", float("nan")) if isinstance(result["sideband_diag_true"], dict) else float("nan"),
                    "twin_leak_decoy": result["sideband_diag_decoy"].get("twin_leak", float("nan")) if isinstance(result["sideband_diag_decoy"], dict) else float("nan"),
                    "winding_true": result.get("winding_true", float("nan")),
                    "winding_decoy": result.get("winding_decoy", float("nan")),
                    "vortex_cx_true": result.get("vortex_cx_true", float("nan")),
                    "vortex_cy_true": result.get("vortex_cy_true", float("nan")),
                    "vortex_cx_decoy": result.get("vortex_cx_decoy", float("nan")),
                    "vortex_cy_decoy": result.get("vortex_cy_decoy", float("nan")),
                    "charge_sum_true": result.get("charge_diag_true", {}).get("sum_charge", float("nan")),
                    "charge_sum_abs_true": result.get("charge_diag_true", {}).get("sum_charge_abs", float("nan")),
                    "charge_n_vort_true": result.get("charge_diag_true", {}).get("n_vort", float("nan")),
                    "charge_cx_true": result.get("charge_diag_true", {}).get("charge_cx", float("nan")),
                    "charge_cy_true": result.get("charge_diag_true", {}).get("charge_cy", float("nan")),
                    "q_float_mean_true": result.get("charge_diag_true", {}).get("q_float_mean", float("nan")),
                    "q_float_std_true": result.get("charge_diag_true", {}).get("q_float_std", float("nan")),
                    "q_float_max_abs_true": result.get("charge_diag_true", {}).get("q_float_max_abs", float("nan")),
                    "n_q_gt_0p25_true": result.get("charge_diag_true", {}).get("n_q_gt_0p25", float("nan")),
                    "n_q_gt_0p5_true": result.get("charge_diag_true", {}).get("n_q_gt_0p5", float("nan")),
                    "q_abs_p50_true": result.get("charge_diag_true", {}).get("q_abs_p50", float("nan")),
                    "q_abs_p90_true": result.get("charge_diag_true", {}).get("q_abs_p90", float("nan")),
                    "n_valid_cells_true": result.get("charge_diag_true", {}).get("n_valid_cells", float("nan")),
                    "valid_cell_frac_true": result.get("charge_diag_true", {}).get("valid_cell_frac", float("nan")),
                    "charge_phase_std_true": result.get("charge_diag_true", {}).get("charge_phase_std", float("nan")),
                    "charge_imag_max_true": result.get("charge_diag_true", {}).get("charge_imag_max", float("nan")),
                    "charge_phase_std_raw_true": result.get("charge_diag_true", {}).get("charge_phase_std_raw", float("nan")),
                    "charge_imag_max_raw_true": result.get("charge_diag_true", {}).get("charge_imag_max_raw", float("nan")),
                    "charge_sum_decoy": result.get("charge_diag_decoy", {}).get("sum_charge", float("nan")),
                    "charge_sum_abs_decoy": result.get("charge_diag_decoy", {}).get("sum_charge_abs", float("nan")),
                    "charge_n_vort_decoy": result.get("charge_diag_decoy", {}).get("n_vort", float("nan")),
                    "charge_cx_decoy": result.get("charge_diag_decoy", {}).get("charge_cx", float("nan")),
                    "charge_cy_decoy": result.get("charge_diag_decoy", {}).get("charge_cy", float("nan")),
                    "q_float_mean_decoy": result.get("charge_diag_decoy", {}).get("q_float_mean", float("nan")),
                    "q_float_std_decoy": result.get("charge_diag_decoy", {}).get("q_float_std", float("nan")),
                    "q_float_max_abs_decoy": result.get("charge_diag_decoy", {}).get("q_float_max_abs", float("nan")),
                    "n_q_gt_0p25_decoy": result.get("charge_diag_decoy", {}).get("n_q_gt_0p25", float("nan")),
                    "n_q_gt_0p5_decoy": result.get("charge_diag_decoy", {}).get("n_q_gt_0p5", float("nan")),
                    "q_abs_p50_decoy": result.get("charge_diag_decoy", {}).get("q_abs_p50", float("nan")),
                    "q_abs_p90_decoy": result.get("charge_diag_decoy", {}).get("q_abs_p90", float("nan")),
                    "n_valid_cells_decoy": result.get("charge_diag_decoy", {}).get("n_valid_cells", float("nan")),
                    "valid_cell_frac_decoy": result.get("charge_diag_decoy", {}).get("valid_cell_frac", float("nan")),
                    "charge_phase_std_decoy": result.get("charge_diag_decoy", {}).get("charge_phase_std", float("nan")),
                    "charge_imag_max_decoy": result.get("charge_diag_decoy", {}).get("charge_imag_max", float("nan")),
                    "charge_phase_std_raw_decoy": result.get("charge_diag_decoy", {}).get("charge_phase_std_raw", float("nan")),
                    "charge_imag_max_raw_decoy": result.get("charge_diag_decoy", {}).get("charge_imag_max_raw", float("nan")),
                    "circ_m_med_true": result.get("circ_true", {}).get("m_circ_med", float("nan")),
                    "circ_m_std_true": result.get("circ_true", {}).get("m_circ_std", float("nan")),
                    "circ_m_n_true": result.get("circ_true", {}).get("m_circ_n", float("nan")),
                    "m_circ_kept_frac_true": result.get("circ_true", {}).get("m_circ_kept_frac", float("nan")),
                    "circ_m_med_decoy": result.get("circ_decoy", {}).get("m_circ_med", float("nan")),
                    "circ_m_std_decoy": result.get("circ_decoy", {}).get("m_circ_std", float("nan")),
                    "circ_m_n_decoy": result.get("circ_decoy", {}).get("m_circ_n", float("nan")),
                    "m_circ_kept_frac_decoy": result.get("circ_decoy", {}).get("m_circ_kept_frac", float("nan")),
                    "proj_ell_-3_true": result.get("proj_true", {}).get("proj_ell_-3", float("nan")),
                    "proj_ell_-2_true": result.get("proj_true", {}).get("proj_ell_-2", float("nan")),
                    "proj_ell_-1_true": result.get("proj_true", {}).get("proj_ell_-1", float("nan")),
                    "proj_ell_1_true": result.get("proj_true", {}).get("proj_ell_1", float("nan")),
                    "proj_ell_2_true": result.get("proj_true", {}).get("proj_ell_2", float("nan")),
                    "proj_ell_3_true": result.get("proj_true", {}).get("proj_ell_3", float("nan")),
                    "proj_ell_-3_decoy": result.get("proj_decoy", {}).get("proj_ell_-3", float("nan")),
                    "proj_ell_-2_decoy": result.get("proj_decoy", {}).get("proj_ell_-2", float("nan")),
                    "proj_ell_-1_decoy": result.get("proj_decoy", {}).get("proj_ell_-1", float("nan")),
                    "proj_ell_1_decoy": result.get("proj_decoy", {}).get("proj_ell_1", float("nan")),
                    "proj_ell_2_decoy": result.get("proj_decoy", {}).get("proj_ell_2", float("nan")),
                    "proj_ell_3_decoy": result.get("proj_decoy", {}).get("proj_ell_3", float("nan")),
                    "proj_peak_true": result.get("proj_peak_true", float("nan")),
                    "proj_peak_ratio_true": result.get("proj_peak_ratio_true", float("nan")),
                    "proj_peak_decoy": result.get("proj_peak_decoy", float("nan")),
                    "proj_peak_ratio_decoy": result.get("proj_peak_ratio_decoy", float("nan")),
                    "proj_center_dx_true": result.get("proj_true", {}).get("proj_center_dx", float("nan")),
                    "proj_center_dy_true": result.get("proj_true", {}).get("proj_center_dy", float("nan")),
                    "proj_center_dx_decoy": result.get("proj_decoy", {}).get("proj_center_dx", float("nan")),
                    "proj_center_dy_decoy": result.get("proj_decoy", {}).get("proj_center_dy", float("nan")),
                    "spec_peak_ell_true": result.get("spec_peak_ell_true", float("nan")),
                    "spec_peak_ratio_true": result.get("spec_peak_ratio_true", float("nan")),
                    "spec_stability_true": result.get("spec_stability_true", float("nan")),
                    "spec_peak_ell_decoy": result.get("spec_peak_ell_decoy", float("nan")),
                    "spec_peak_ratio_decoy": result.get("spec_peak_ratio_decoy", float("nan")),
                    "spec_stability_decoy": result.get("spec_stability_decoy", float("nan")),
                    "spec_med_ell_-3_true": result.get("spec_med_ell_-3_true", float("nan")),
                    "spec_med_ell_-2_true": result.get("spec_med_ell_-2_true", float("nan")),
                    "spec_med_ell_-1_true": result.get("spec_med_ell_-1_true", float("nan")),
                    "spec_med_ell_1_true": result.get("spec_med_ell_1_true", float("nan")),
                    "spec_med_ell_2_true": result.get("spec_med_ell_2_true", float("nan")),
                    "spec_med_ell_3_true": result.get("spec_med_ell_3_true", float("nan")),
                    "spec_std_ell_-3_true": result.get("spec_std_ell_-3_true", float("nan")),
                    "spec_std_ell_-2_true": result.get("spec_std_ell_-2_true", float("nan")),
                    "spec_std_ell_-1_true": result.get("spec_std_ell_-1_true", float("nan")),
                    "spec_std_ell_1_true": result.get("spec_std_ell_1_true", float("nan")),
                    "spec_std_ell_2_true": result.get("spec_std_ell_2_true", float("nan")),
                    "spec_std_ell_3_true": result.get("spec_std_ell_3_true", float("nan")),
                    "spec_med_ell_-3_decoy": result.get("spec_med_ell_-3_decoy", float("nan")),
                    "spec_med_ell_-2_decoy": result.get("spec_med_ell_-2_decoy", float("nan")),
                    "spec_med_ell_-1_decoy": result.get("spec_med_ell_-1_decoy", float("nan")),
                    "spec_med_ell_1_decoy": result.get("spec_med_ell_1_decoy", float("nan")),
                    "spec_med_ell_2_decoy": result.get("spec_med_ell_2_decoy", float("nan")),
                    "spec_med_ell_3_decoy": result.get("spec_med_ell_3_decoy", float("nan")),
                    "spec_std_ell_-3_decoy": result.get("spec_std_ell_-3_decoy", float("nan")),
                    "spec_std_ell_-2_decoy": result.get("spec_std_ell_-2_decoy", float("nan")),
                    "spec_std_ell_-1_decoy": result.get("spec_std_ell_-1_decoy", float("nan")),
                    "spec_std_ell_1_decoy": result.get("spec_std_ell_1_decoy", float("nan")),
                    "spec_std_ell_2_decoy": result.get("spec_std_ell_2_decoy", float("nan")),
                    "spec_std_ell_3_decoy": result.get("spec_std_ell_3_decoy", float("nan")),
                    "oam_peak_pre_true": result.get("oam_pre_true", {}).get("m_peak", float("nan")),
                    "oam_peak_ratio_pre_true": result.get("oam_pre_true", {}).get("m_peak_ratio", float("nan")),
                    "oam_pre_raw_peak_true": result.get("oam_pre_true", {}).get("m_peak_raw", float("nan")),
                    "oam_pre_raw_ratio_true": result.get("oam_pre_true", {}).get("m_peak_ratio_raw", float("nan")),
                    "oam_pre_gated_true": int(bool(result.get("oam_pre_true", {}).get("gated", False))),
                    "oam_pre_gate_reason_true": result.get("oam_pre_true", {}).get("gate_reason", ""),
                    "oam_peak_pre_decoy": result.get("oam_pre_decoy", {}).get("m_peak", float("nan")),
                    "oam_peak_ratio_pre_decoy": result.get("oam_pre_decoy", {}).get("m_peak_ratio", float("nan")),
                    "oam_pre_raw_peak_decoy": result.get("oam_pre_decoy", {}).get("m_peak_raw", float("nan")),
                    "oam_pre_raw_ratio_decoy": result.get("oam_pre_decoy", {}).get("m_peak_ratio_raw", float("nan")),
                    "oam_pre_gated_decoy": int(bool(result.get("oam_pre_decoy", {}).get("gated", False))),
                    "oam_pre_gate_reason_decoy": result.get("oam_pre_decoy", {}).get("gate_reason", ""),
                    "oam_peak_post_true": result.get("oam_peak_post_true", float("nan")),
                    "oam_peak_ratio_post_true": result.get("oam_peak_ratio_post_true", float("nan")),
                    "oam_peak_post_decoy": result.get("oam_peak_post_decoy", float("nan")),
                    "oam_peak_ratio_post_decoy": result.get("oam_peak_ratio_post_decoy", float("nan")),
                    "knee_amp_true": result.get("knee_amp_true", float("nan")),
                    "knee_amp_decoy": result.get("knee_amp_decoy", float("nan")),
                    "knee_phasegrad_true": result.get("knee_phasegrad_true", float("nan")),
                    "knee_phasegrad_decoy": result.get("knee_phasegrad_decoy", float("nan")),
                    "leak_penalty": result.get("leak_penalty", float("nan")),
                    "leak_fail": int(result.get("leak_fail", False)),
                    "peak_y_true": result.get("peak_y_true", float("nan")),
                    "peak_x_true": result.get("peak_x_true", float("nan")),
                    "crop_radius_true": result.get("crop_radius_true", float("nan")),
                    "peak_y_decoy": result.get("peak_y_decoy", float("nan")),
                    "peak_x_decoy": result.get("peak_x_decoy", float("nan")),
                    "crop_radius_decoy": result.get("crop_radius_decoy", float("nan")),
                    "peak_y_pred_true": result["sideband_diag_true"].get("peak_y_pred", float("nan")) if isinstance(result["sideband_diag_true"], dict) else float("nan"),
                    "peak_x_pred_true": result["sideband_diag_true"].get("peak_x_pred", float("nan")) if isinstance(result["sideband_diag_true"], dict) else float("nan"),
                    "peak_y_pred_decoy": result["sideband_diag_decoy"].get("peak_y_pred", float("nan")) if isinstance(result["sideband_diag_decoy"], dict) else float("nan"),
                    "peak_x_pred_decoy": result["sideband_diag_decoy"].get("peak_x_pred", float("nan")) if isinstance(result["sideband_diag_decoy"], dict) else float("nan"),
                    "r_min_true": result["sideband_diag_true"].get("r_min", float("nan")) if isinstance(result["sideband_diag_true"], dict) else float("nan"),
                    "r_min_decoy": result["sideband_diag_decoy"].get("r_min", float("nan")) if isinstance(result["sideband_diag_decoy"], dict) else float("nan"),
                    "dist_detect_pred_true": result["sideband_diag_true"].get("dist_detect_pred", float("nan")) if isinstance(result["sideband_diag_true"], dict) else float("nan"),
                    "dist_detect_pred_decoy": result["sideband_diag_decoy"].get("dist_detect_pred", float("nan")) if isinstance(result["sideband_diag_decoy"], dict) else float("nan"),
                    "kx_used_bins_true": result["sideband_diag_true"].get("kx_used_bins", float("nan")) if isinstance(result["sideband_diag_true"], dict) else float("nan"),
                    "ky_used_bins_true": result["sideband_diag_true"].get("ky_used_bins", float("nan")) if isinstance(result["sideband_diag_true"], dict) else float("nan"),
                    "kx_used_bins_decoy": result["sideband_diag_decoy"].get("kx_used_bins", float("nan")) if isinstance(result["sideband_diag_decoy"], dict) else float("nan"),
                    "ky_used_bins_decoy": result["sideband_diag_decoy"].get("ky_used_bins", float("nan")) if isinstance(result["sideband_diag_decoy"], dict) else float("nan"),
                    "dist_detect_pred_plus_true": result["sideband_diag_true"].get("dist_detect_pred_plus", float("nan")) if isinstance(result["sideband_diag_true"], dict) else float("nan"),
                    "dist_detect_pred_minus_true": result["sideband_diag_true"].get("dist_detect_pred_minus", float("nan")) if isinstance(result["sideband_diag_true"], dict) else float("nan"),
                    "dist_detect_pred_plus_decoy": result["sideband_diag_decoy"].get("dist_detect_pred_plus", float("nan")) if isinstance(result["sideband_diag_decoy"], dict) else float("nan"),
                    "dist_detect_pred_minus_decoy": result["sideband_diag_decoy"].get("dist_detect_pred_minus", float("nan")) if isinstance(result["sideband_diag_decoy"], dict) else float("nan"),
                }
                for key, val in result["metrics_true"].items():
                    row[f"true_{key}"] = val
                for key, val in result["metrics_decoy"].items():
                    row[f"decoy_{key}"] = val
                rows.append(row)

                if save_images:
                    stem = images_dir / f"trial{trial}_{name}_{approach_lower}"
                    cv2.imwrite(str(stem.with_name(f"{stem.name}_true_scrambled.png")), (result["true_scrambled"] * 255).astype(np.uint8))
                    cv2.imwrite(str(stem.with_name(f"{stem.name}_decoy_scrambled.png")), (result["decoy_scrambled"] * 255).astype(np.uint8))
                    cv2.imwrite(str(stem.with_name(f"{stem.name}_decoded_true.png")), (result["decoded_true"] * 255).astype(np.uint8))
                    cv2.imwrite(str(stem.with_name(f"{stem.name}_decoded_decoy.png")), (result["decoded_decoy"] * 255).astype(np.uint8))

    df = pd.DataFrame(rows)
    for path in (csv_path, stamped):
        try:
            df.to_csv(path, index=False)
        except PermissionError:
            continue
    print(f"Saved results to {csv_path} and {stamped}")
    return df


# --------------------------
# CLI
# --------------------------
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Spiral/OAM-informed decoding experiment with randomized controls."
    )
    parser.add_argument("--out-dir", type=Path, required=True, help="Output directory for CSV and optional PNGs.")
    parser.add_argument("--trials", type=int, default=3, help="Number of repetitions per pattern.")
    parser.add_argument("--image-size", type=int, default=256, help="Generated pattern size (pixels).")
    parser.add_argument("--gaussian-sigma", type=float, default=0.02, help="Gaussian noise sigma added after encoding.")
    parser.add_argument("--poisson-scale", type=float, default=0.0, help="Optional Poisson noise scaling (0 to disable).")
    parser.add_argument(
        "--save-images", action=argparse.BooleanOptionalAction, default=True, help="Save scrambled/decoded PNGs."
    )
    parser.add_argument("--seed", type=int, default=7, help="RNG seed for reproducibility.")
    parser.add_argument(
        "--patterns",
        type=str,
        default="ring,bars,spiral_hint,comma",
        help="Comma list of patterns to include (subset of generate_patterns keys).",
    )
    parser.add_argument("--ref-amp", type=float, default=0.0, help="Reference wave amplitude for interferometric measurement.")
    parser.add_argument("--ref-kx", type=float, default=0.0, help="Reference wave tilt kx.")
    parser.add_argument("--ref-ky", type=float, default=0.0, help="Reference wave tilt ky.")
    parser.add_argument(
        "--carrier-radius",
        type=float,
        default=0.0,
        help="Optional carrier radius (cycles/image) to derive ref_kx/ref_ky using carrier-ratio.",
    )
    parser.add_argument(
        "--carrier-ratio",
        type=float,
        default=np.pi / np.e,
        help="kx/ky ratio for carrier direction (default pi/e).",
    )
    parser.add_argument(
        "--ells",
        type=str,
        default="-3,-2,-1,1,2,3",
        help="Comma list of ell values for candidate grid (0 is ignored).",
    )
    parser.add_argument(
        "--pitches",
        type=str,
        default="0.0,1.0,1.5,2.5,3.0,3.5",
        help="Comma list of pitch values for candidate grid.",
    )
    parser.add_argument(
        "--knees",
        type=str,
        default="0.0,0.5,1.0,1.5",
        help="Comma list of knee values for candidate grid.",
    )
    parser.add_argument(
        "--approaches",
        type=str,
        default="baseline,interferometric,angular",
        help="Comma list of decoding approaches to run (baseline, interferometric, angular).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    noise = NoiseParams(gaussian_sigma=args.gaussian_sigma, poisson_scale=args.poisson_scale)
    pattern_filter = {p.strip() for p in args.patterns.split(",") if p.strip()}
    approaches = [a.strip() for a in args.approaches.split(",") if a.strip()]
    def parse_list(raw: str, cast: Callable) -> List:
        vals = []
        for part in raw.split(","):
            part = part.strip()
            if not part:
                continue
            try:
                vals.append(cast(part))
            except Exception:
                continue
        return vals
    grid_ells = parse_list(args.ells, int)
    grid_pitches = parse_list(args.pitches, float)
    grid_knees = parse_list(args.knees, float)
    ref_kx = args.ref_kx
    ref_ky = args.ref_ky
    if args.carrier_radius > 0:
        kx_c, ky_c = carrier_from_ratio(args.carrier_radius, args.carrier_ratio)
        ref_kx = kx_c
        ref_ky = ky_c
    df = run_experiment(
        out_dir=args.out_dir,
        trials=args.trials,
        image_size=args.image_size,
        noise=noise,
        save_images=args.save_images,
        seed=args.seed,
        pattern_filter=pattern_filter,
        ref_amp=args.ref_amp,
        ref_kx=ref_kx,
        ref_ky=ref_ky,
        approaches=approaches,
        ells=grid_ells,
        pitches=grid_pitches,
        knees=grid_knees,
    )
    summary = df.describe().transpose()[["mean", "std"]]
    print("Saved results to", args.out_dir / "spiral_decoder_results.csv")
    print("Metric summary (mean/std):")
    print(summary[summary.index.str.contains("ssim|ncc|psnr")])


if __name__ == "__main__":
    main()
