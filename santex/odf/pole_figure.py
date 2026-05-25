"""
Pole figure and inverse pole figure utilities.

Pole figure  (PF):  projects a crystal direction h into the sample frame.
Inverse PF   (IPF): projects a sample direction r into the crystal frame.
"""

from __future__ import annotations
import numpy as np


# ---------------------------------------------------------------------------
# Stereographic projection helpers
# ---------------------------------------------------------------------------

def stereo_project(vectors: np.ndarray, hemisphere: str = "upper") -> np.ndarray:
    """
    Stereographic projection from the unit sphere to the equatorial plane.

    Parameters
    ----------
    vectors   : (N, 3) unit vectors
    hemisphere: 'upper' (z ≥ 0) or 'lower'; vectors in the wrong hemisphere
                are negated (antipodal mapping).

    Returns
    -------
    xy : (N, 2) projected coordinates
    """
    v = np.asarray(vectors, dtype=float).copy()
    if hemisphere == "upper":
        neg = v[:, 2] < 0
        v[neg] = -v[neg]
    else:
        pos = v[:, 2] > 0
        v[pos] = -v[pos]

    denom = 1.0 + v[:, 2]
    # Guard division by zero (south pole → infinity)
    denom = np.where(np.abs(denom) < 1e-12, 1e-12, denom)
    X = v[:, 0] / denom
    Y = v[:, 1] / denom
    return np.stack([X, Y], axis=1)


def _rotation_matrices_from_euler(phi1, PHI, phi2):
    """Vectorised Bunge ZXZ → (N, 3, 3)."""
    from santex.grains.reconstruction import euler_to_rotation_matrix
    phi1 = np.asarray(phi1, dtype=float)
    PHI  = np.asarray(PHI,  dtype=float)
    phi2 = np.asarray(phi2, dtype=float)
    return euler_to_rotation_matrix(phi1, PHI, phi2)


# ---------------------------------------------------------------------------
# Crystal symmetry helpers
# ---------------------------------------------------------------------------

def _symmetry_ops(crystal_symmetry: str) -> np.ndarray:
    from santex.grains.symmetry import SYMMETRY
    return SYMMETRY.get(crystal_symmetry.lower(), SYMMETRY["triclinic"])


def _sample_symmetry_ops(sample_symmetry: str) -> np.ndarray:
    """Sample symmetry operators (orthorhombic is the most common in geology)."""
    if sample_symmetry.lower() in ("orthorhombic", "mmm", "d2h"):
        # 4 proper rotations: identity, 180° around X, Y, Z
        return np.array([
            np.eye(3),
            np.diag([ 1, -1, -1]),  # 180° Z
            np.diag([-1,  1, -1]),  # 180° Y
            np.diag([-1, -1,  1]),  # 180° X
        ], dtype=float)
    # Default: triclinic (identity only)
    return np.array([np.eye(3)], dtype=float)


# ---------------------------------------------------------------------------
# Pole figure: project crystal direction into sample frame
# ---------------------------------------------------------------------------

def calc_pole_figure(euler_deg: np.ndarray,
                     hkl: tuple[float, float, float] = (0, 0, 1),
                     crystal_symmetry: str = "cubic",
                     sample_symmetry: str = "none",
                     weights: np.ndarray | None = None
                     ) -> tuple[np.ndarray, np.ndarray | None]:
    """
    Compute pole figure vectors for a given crystal direction {hkl}.

    For each orientation g and each crystal-symmetry operator s_c:
        r = g @ s_c @ h   →  unit vector in sample frame

    Additionally, for each sample-symmetry operator s_s: r' = s_s @ r.

    Parameters
    ----------
    euler_deg         : (N, 3) Euler angles in degrees
    hkl               : crystal direction (need not be normalised)
    crystal_symmetry  : 'cubic', 'hexagonal', 'orthorhombic', 'triclinic', …
    sample_symmetry   : 'none', 'orthorhombic', …
    weights           : (N,) per-orientation weights

    Returns
    -------
    vectors : (M, 3) unit vectors in sample frame (M = N × N_csym × N_ssym)
    w_out   : (M,) weights or None
    """
    euler_rad = np.radians(euler_deg)
    R = _rotation_matrices_from_euler(
        euler_rad[:, 0], euler_rad[:, 1], euler_rad[:, 2]
    )  # (N, 3, 3)

    h = np.array(hkl, dtype=float)
    h = h / (np.linalg.norm(h) + 1e-15)

    csym = _symmetry_ops(crystal_symmetry)   # (S_c, 3, 3)
    ssym = _sample_symmetry_ops(sample_symmetry)  # (S_s, 3, 3)

    # Apply crystal symmetry to h: (S_c, 3)
    h_sym = (csym @ h)  # (S_c, 3)

    # Rotate each h_sym by each orientation: r[n, sc] = R[n] @ h_sym[sc]
    # (N, S_c, 3)
    r = np.einsum("nij,sj->nsi", R, h_sym)

    # Apply sample symmetry: r_full[n, sc, ss] = ssym[ss] @ r[n, sc]
    if len(ssym) > 1:
        r = np.einsum("sij,nkj->nksi", ssym, r)   # (N, S_c, S_s, 3)
        r = r.reshape(-1, 3)
    else:
        r = r.reshape(-1, 3)

    # Normalize
    norms = np.linalg.norm(r, axis=1, keepdims=True)
    r = r / np.where(norms < 1e-15, 1.0, norms)

    # Expand weights
    w_out = None
    if weights is not None:
        n_repeat = len(r) // len(weights)
        w_out = np.repeat(weights, n_repeat)

    return r, w_out


# ---------------------------------------------------------------------------
# Inverse pole figure: project sample direction into crystal frame
# ---------------------------------------------------------------------------

def calc_inverse_pole_figure(euler_deg: np.ndarray,
                              sample_dir: tuple[float, float, float] = (0, 0, 1),
                              crystal_symmetry: str = "cubic"
                              ) -> np.ndarray:
    """
    Compute crystal-frame directions corresponding to a given sample direction.

    For each orientation g:  c = g^T @ r  (crystal frame direction)
    Then apply crystal symmetry to bring c into the fundamental sector.

    Parameters
    ----------
    euler_deg        : (N, 3) Euler angles in degrees
    sample_dir       : unit vector in sample frame (e.g. [0,0,1] = ND)
    crystal_symmetry : symmetry string

    Returns
    -------
    crystal_dirs : (N, 3) unit vectors in crystal frame, in fundamental sector
    """
    euler_rad = np.radians(euler_deg)
    R = _rotation_matrices_from_euler(
        euler_rad[:, 0], euler_rad[:, 1], euler_rad[:, 2]
    )  # (N, 3, 3)

    r = np.array(sample_dir, dtype=float)
    r = r / (np.linalg.norm(r) + 1e-15)

    # Crystal direction: c = R^T @ r = R^{-1} @ r
    c = np.einsum("nji,j->ni", R, r)  # (N, 3)

    # Bring into fundamental sector using crystal symmetry
    csym = _symmetry_ops(crystal_symmetry)
    c = _into_fundamental_sector(c, csym, crystal_symmetry)

    return c


def _into_fundamental_sector(c: np.ndarray, csym: np.ndarray,
                              crystal_symmetry: str) -> np.ndarray:
    """
    Apply crystal symmetry operators to bring vectors into the fundamental sector.

    For cubic (m-3m): fundamental sector is h ≥ k ≥ l ≥ 0.
    For other symmetries: use first positive octant approximation.
    """
    s = crystal_symmetry.lower()

    if "cubic" in s or "m-3m" in s or "m3m" in s:
        # For cubic: take absolute values and sort descending (h ≥ k ≥ l ≥ 0)
        c_abs = np.abs(c)
        c_abs = np.sort(c_abs, axis=1)[:, ::-1]   # descending
        norms = np.linalg.norm(c_abs, axis=1, keepdims=True)
        return c_abs / np.where(norms < 1e-15, 1.0, norms)

    # General: apply all sym ops, pick the one with max z (closest to north pole)
    best = c.copy()
    for s_op in csym:
        c_s = (s_op @ c.T).T  # (N, 3)
        # Prefer vectors in positive hemisphere with max z
        mask = c_s[:, 2] > best[:, 2]
        best[mask] = c_s[mask]
    norms = np.linalg.norm(best, axis=1, keepdims=True)
    return best / np.where(norms < 1e-15, 1.0, norms)


# ---------------------------------------------------------------------------
# IPF colour coding
# ---------------------------------------------------------------------------

def ipf_color(crystal_dirs: np.ndarray,
              crystal_symmetry: str = "cubic") -> np.ndarray:
    """
    Map crystal-frame unit vectors (in fundamental sector) to RGB colours.

    Colour assignment (cubic):
      [001] = Blue,  [011] = Green,  [111] = Red

    Parameters
    ----------
    crystal_dirs  : (N, 3) unit vectors already in fundamental sector
    crystal_symmetry : for lookup of colour key (currently cubic implemented)

    Returns
    -------
    rgb : (N, 3) array of RGB values in [0, 1]
    """
    d = np.asarray(crystal_dirs, dtype=float)
    norms = np.linalg.norm(d, axis=1, keepdims=True)
    d = d / np.where(norms < 1e-15, 1.0, norms)

    s = crystal_symmetry.lower()

    if "cubic" in s or "m-3m" in s or "m3m" in s:
        return _ipf_color_cubic(d)
    if "hex" in s or "6/mmm" in s:
        return _ipf_color_hex(d)
    # Generic: map z→B, y→G, x→R
    rgb = np.abs(d)
    mx = rgb.max(axis=1, keepdims=True)
    return rgb / np.where(mx < 1e-15, 1.0, mx)


def _ipf_color_cubic(d: np.ndarray) -> np.ndarray:
    """
    Cubic IPF coloring.
    d: (N,3) in fundamental domain, h >= k >= l >= 0 (sorted descending).

    Corners:
      [001]=(0,0,1) → Blue
      [011]=(0,1,1)/√2 → Green
      [111]=(1,1,1)/√3 → Red

    Linear interpolation:
      R = h * √3              (max 1 at [111])
      G = (k - h) * √2       (max 1 at [011] where h=0, k=l=1/√2)
      B = l - k               (max 1 at [001] where k=0, l=1)
    """
    h, k, l = d[:, 0], d[:, 1], d[:, 2]

    R = np.clip(h * np.sqrt(3), 0, 1)
    G = np.clip((k - h) * np.sqrt(2), 0, 1)
    B = np.clip(l - k, 0, 1)

    # Normalise so the brightest channel = 1
    rgb = np.stack([R, G, B], axis=1)
    mx = rgb.max(axis=1, keepdims=True)
    rgb = rgb / np.where(mx < 1e-15, 1.0, mx)
    return np.clip(rgb, 0, 1)


def _ipf_color_hex(d: np.ndarray) -> np.ndarray:
    """
    Hexagonal IPF coloring (6/mmm).
    Corners: [0001]→Blue, [10-10]→Red, [2-1-10]→Green.
    d assumed to be in fundamental sector (z ≥ 0, x ≥ 0, in first sextant).
    """
    # Map (x, y, z) — approximate coloring
    R = np.clip(d[:, 0], 0, 1)  # [10-10] character
    G = np.clip(d[:, 1], 0, 1)  # [01-10] character
    B = np.clip(d[:, 2], 0, 1)  # [0001] character
    rgb = np.stack([R, G, B], axis=1)
    mx = rgb.max(axis=1, keepdims=True)
    return rgb / np.where(mx < 1e-15, 1.0, mx)


# ---------------------------------------------------------------------------
# IPF map: colour each EBSD pixel
# ---------------------------------------------------------------------------

def ipf_map_colors(euler_deg: np.ndarray,
                   sample_dir: tuple[float, float, float] = (0, 0, 1),
                   crystal_symmetry: str = "cubic") -> np.ndarray:
    """
    Compute per-pixel IPF colours for an EBSD dataset.

    Parameters
    ----------
    euler_deg        : (N, 3) Euler angles in degrees
    sample_dir       : reference sample direction (default: ND = [0,0,1])
    crystal_symmetry : crystal symmetry string

    Returns
    -------
    rgb : (N, 3) RGB values in [0, 1]
    """
    crystal_dirs = calc_inverse_pole_figure(euler_deg, sample_dir, crystal_symmetry)
    return ipf_color(crystal_dirs, crystal_symmetry)


# ---------------------------------------------------------------------------
# Sphere KDE for smooth pole figure contours
# ---------------------------------------------------------------------------

def sphere_kde(vectors: np.ndarray, eval_pts: np.ndarray,
               bandwidth: float = 0.1,
               weights: np.ndarray | None = None) -> np.ndarray:
    """
    Kernel density estimate on the unit sphere (S²).

    Uses a von Mises–Fisher kernel: K(v·u) = exp(κ(v·u - 1))
    where κ = 1 / bandwidth².

    Parameters
    ----------
    vectors   : (N, 3) data unit vectors
    eval_pts  : (M, 3) evaluation unit vectors
    bandwidth : angular bandwidth in radians
    weights   : (N,) optional weights

    Returns
    -------
    density : (M,) density values
    """
    kappa = 1.0 / (bandwidth ** 2 + 1e-15)
    # dot products: (M, N)
    dots = np.clip(eval_pts @ vectors.T, -1.0, 1.0)
    K = np.exp(kappa * (dots - 1.0))   # (M, N)
    if weights is not None:
        K = K * weights[None, :]
        density = K.sum(axis=1) / (weights.sum() + 1e-15)
    else:
        density = K.mean(axis=1)
    return density


def pole_figure_grid(n_pts: int = 500) -> tuple[np.ndarray, np.ndarray]:
    """
    Generate a regular grid on the upper hemisphere for PF contour plotting.

    Returns
    -------
    xy     : (M, 2) stereographic coordinates
    vecs   : (M, 3) corresponding unit vectors
    """
    # Fibonacci spiral on sphere
    golden = (1 + np.sqrt(5)) / 2
    i = np.arange(n_pts)
    theta = 2 * np.pi * i / golden        # azimuth
    phi   = np.arccos(1 - 2 * (i + 0.5) / n_pts)  # colatitude

    x = np.sin(phi) * np.cos(theta)
    y = np.sin(phi) * np.sin(theta)
    z = np.cos(phi)
    vecs = np.stack([x, y, z], axis=1)

    # Keep upper hemisphere
    mask = vecs[:, 2] >= 0
    vecs = vecs[mask]
    xy = stereo_project(vecs)
    return xy, vecs
