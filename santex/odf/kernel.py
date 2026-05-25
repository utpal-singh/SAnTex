"""
ODF kernel functions for kernel density estimation on SO(3).

Reference: Hielscher & Schaeben (2008), J. Appl. Cryst. 41, 1024-1037.
"""

from __future__ import annotations
import numpy as np


# ---------------------------------------------------------------------------
# de la Vallée-Poussin kernel (most common in EBSD/texture analysis)
# ---------------------------------------------------------------------------

def de_la_vallee_poussin(omega: np.ndarray, kappa: float) -> np.ndarray:
    """
    de la Vallée-Poussin kernel evaluated at misorientation angles.

    K_κ(ω) = C_κ · cos^(2κ)(ω/2)

    Normalization: C_κ = π / (2 · B(κ + ½, 3/2))  where B is the Beta function.
    This ensures ∫_{SO(3)} K_κ(ω(g, g₀)) dg = 1, so a random ODF gives f ≈ 1 m.r.d.

    Reference: Hielscher & Schaeben (2008), eq. (5).

    Parameters
    ----------
    omega : misorientation angles in radians, shape (...)
    kappa : sharpness parameter (larger = sharper; ≈91 for 10° HWHM)

    Returns
    -------
    values : kernel values, same shape as omega
    """
    from scipy.special import betaln
    kappa = float(kappa)
    cos_half = np.cos(np.asarray(omega, dtype=np.float64) / 2.0)
    # C_κ = π / (2 B(κ+½, 3/2)); use log-space to avoid overflow for large κ
    log_C = np.log(np.pi) - np.log(2.0) - betaln(kappa + 0.5, 1.5)
    C = np.exp(log_C)
    return C * np.maximum(cos_half, 0.0) ** (2.0 * kappa)


def von_mises_fisher(omega: np.ndarray, kappa: float) -> np.ndarray:
    """
    von Mises–Fisher kernel on SO(3).

    K_κ(ω) ∝ exp(κ · cos(ω))

    Parameters
    ----------
    omega : misorientation angles in radians
    kappa : concentration parameter (larger = sharper)
    """
    return np.exp(kappa * (np.cos(np.asarray(omega, dtype=np.float64)) - 1.0))


def gaussian_kernel(omega: np.ndarray, halfwidth_rad: float) -> np.ndarray:
    """Wrapped Gaussian kernel (approximate, fast)."""
    sigma2 = halfwidth_rad ** 2 / (8.0 * np.log(2.0))
    return np.exp(-0.5 * np.asarray(omega, dtype=np.float64) ** 2 / sigma2)


# ---------------------------------------------------------------------------
# Bandwidth helpers
# ---------------------------------------------------------------------------

def halfwidth_to_kappa(halfwidth_deg: float) -> float:
    """Convert FWHM half-width (degrees) to de la Vallée-Poussin κ."""
    hw = np.radians(halfwidth_deg)
    cos_val = np.cos(hw / 2.0)
    if cos_val <= 0.0:
        return 1.0
    kappa = np.log(0.5) / (2.0 * np.log(cos_val + 1e-15))
    return max(1.0, kappa)


def kappa_to_halfwidth(kappa: float) -> float:
    """Convert de la Vallée-Poussin κ to FWHM half-width in degrees."""
    cos_val = 0.5 ** (1.0 / (2.0 * kappa))
    return float(np.degrees(2.0 * np.arccos(np.clip(cos_val, -1, 1))))
