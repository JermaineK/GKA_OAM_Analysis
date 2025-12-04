"""
Symmetry-driven filters for OAM / GKA analysis.

These are "projectors onto symmetries" rather than aesthetic tweaks:

Filter A: Parity–odd ring contrast (eta_odd)
    I_even(r,theta) = 0.5*(I + I(theta+pi))
    I_odd (r,theta) = 0.5*(I - I(theta+pi))
    eta_odd(r) = |mean_theta I_odd| / (mean_theta I_even + eps)

Filter B: Spiral harmonic contrast (eta_spiral)
    Keep low angular modes (m=1,2) and their odd (imag) parts.
    P_odd = Im(a1)^2 + Im(a2)^2, where a_m = <I * e^{-im theta}>
    eta_spiral(r) = P_odd / (total_power + eps)

Filter C: Fixed log-spiral kernel projection (eta_kernel)
    K_alpha(r,theta) = cos(alpha * ln r - theta), alpha fixed (theory-driven).
    F_alpha(r) = <I * K_alpha>_theta; eta_kernel = |F_alpha| / (mean_even + eps)

Filter D: Slow-tick projection for time series (eta_slow)
    eta_slow = |Fourier component at f_star|.

All filters are linear in the data and use only fixed parameters (eps, alpha, f_star).
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Tuple

import numpy as np
from scipy.signal import savgol_filter


def parity_odd_contrast_ring(I_r_theta: np.ndarray, eps: float = 1e-9) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute parity-odd contrast per radius.

    Parameters
    ----------
    I_r_theta : array, shape (n_r, n_theta)
        Intensity on a polar grid sampled uniformly in theta.
    eps : float
        Small constant to avoid division by zero.

    Returns
    -------
    eta_odd : array, shape (n_r,)
        Parity-odd contrast per radius.
    mean_even : array, shape (n_r,)
        Mean even component per radius (can be used for normalization).
    """
    n_theta = I_r_theta.shape[1]
    I_roll = np.roll(I_r_theta, shift=n_theta // 2, axis=1)  # theta -> theta+pi
    I_even = 0.5 * (I_r_theta + I_roll)
    I_odd = 0.5 * (I_r_theta - I_roll)
    mean_even = I_even.mean(axis=1)
    mean_odd = I_odd.mean(axis=1)
    eta_odd = np.abs(mean_odd) / (mean_even + eps)
    return eta_odd, mean_even


def spiral_harmonic_contrast(I_r_theta: np.ndarray, eps: float = 1e-9) -> np.ndarray:
    """
    Spiral harmonic contrast using m=1,2 angular Fourier modes (imag parts).

    Returns eta_spiral(r) = (Im a1)^2 + (Im a2)^2 / (total power + eps).
    """
    n_theta = I_r_theta.shape[1]
    theta = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
    basis1 = np.exp(-1j * theta)
    basis2 = np.exp(-2j * theta)
    a1 = (I_r_theta * basis1).mean(axis=1)
    a2 = (I_r_theta * basis2).mean(axis=1)
    p_odd = (a1.imag**2 + a2.imag**2)
    p_tot = (np.abs(I_r_theta) ** 2).mean(axis=1)
    return p_odd / (p_tot + eps)


def log_spiral_kernel(r_vals: np.ndarray, theta_vals: np.ndarray, alpha: float) -> np.ndarray:
    """
    Build a fixed log-spiral kernel K_alpha(r,theta) = cos(alpha ln r - theta).
    """
    R, T = np.meshgrid(r_vals, theta_vals, indexing="ij")
    return np.cos(alpha * np.log(R + 1e-9) - T)


def kernel_projection(I_r_theta: np.ndarray, r_vals: np.ndarray, alpha: float, eps: float = 1e-9) -> Tuple[np.ndarray, np.ndarray]:
    """
    Project onto a fixed log-spiral kernel with pitch alpha.

    Returns
    -------
    eta_kernel : array (n_r,)
        |<I*K>| / (mean_even + eps)
    mean_even : array (n_r,)
        Same as in parity_odd_contrast_ring.
    """
    n_theta = I_r_theta.shape[1]
    theta = np.linspace(0, 2 * np.pi, n_theta, endpoint=False)
    K = log_spiral_kernel(r_vals, theta, alpha)
    proj = (I_r_theta * K).mean(axis=1)
    eta_even, mean_even = parity_odd_contrast_ring(I_r_theta, eps)  # reuse even part
    eta_kernel = np.abs(proj) / (mean_even + eps)
    return eta_kernel, mean_even


def slow_tick_projection(series: np.ndarray, dt: float, f_star: float) -> float:
    """
    Projection of a time series onto a single sinusoid at f_star.

    Returns magnitude of Fourier component at f_star (eta_slow).
    """
    series = np.asarray(series, float)
    n = len(series)
    t = np.arange(n) * dt
    omega = 2 * np.pi * f_star
    A = (2.0 / n) * np.sum(series * np.cos(omega * t))
    B = (2.0 / n) * np.sum(series * np.sin(omega * t))
    return float(np.sqrt(A**2 + B**2))


@dataclass
class KneeResult:
    r_knee: float
    d2: np.ndarray
    log_r: np.ndarray
    log_eta: np.ndarray


def knee_on_loglog(r: np.ndarray, eta: np.ndarray, window: int = 9, polyorder: int = 2) -> KneeResult:
    """
    Knee via curvature on log-log: smooth log(eta) vs log(r) with Savitzky–Golay, then find max |d2|.
    """
    mask = (eta > 0) & (r > 0) & np.isfinite(eta) & np.isfinite(r)
    log_r = np.log(r[mask])
    log_eta = np.log(eta[mask])
    if log_r.size < window:
        window = max(5, 2 * (log_r.size // 2) + 1)
    d1 = np.gradient(log_eta, log_r)
    d1s = savgol_filter(d1, window_length=max(3, window), polyorder=polyorder)
    d2 = np.gradient(d1s, log_r)
    idx = int(np.argmax(np.abs(d2)))
    return KneeResult(r_knee=float(np.exp(log_r[idx])), d2=d2, log_r=log_r, log_eta=log_eta)
