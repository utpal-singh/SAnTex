"""
Pixel-level orientation analysis for EBSD grain maps.

Functions
---------
calc_kam(grains, kernel_size, max_angle_deg)
    Kernel Average Misorientation — mean misorientation of each pixel with
    its immediate neighbours (angle < max_angle_deg).

calc_grod(grains)
    Grain Reference Orientation Deviation — misorientation of each pixel from
    the mean orientation of its grain.

calc_misorientation_distribution(grain_boundary, bins)
    Binned misorientation angle distribution (length-weighted).

calc_odf_texture_index(euler_angles_deg)
    Simple texture index from orientation data (M-index approximation).
"""

from __future__ import annotations
import numpy as np
import pandas as pd


def calc_kam(grains, kernel_size: int = 1,
             max_angle_deg: float = 5.0) -> np.ndarray:
    """
    Kernel Average Misorientation (KAM) per pixel.

    Parameters
    ----------
    grains        : Grain2D object
    kernel_size   : number of pixel steps to consider as neighbours (default 1 = 4-conn)
    max_angle_deg : maximum misorientation to include in the average

    Returns
    -------
    kam : (N,) float array (degrees), NaN for unindexed pixels
    """
    from .reconstruction import euler_to_rotation_matrix, build_pixel_adjacency

    data = grains.data
    g    = grains.grain_ids
    n    = len(data)

    phi1 = np.radians(data['Euler1'].values)
    PHI  = np.radians(data['Euler2'].values)
    phi2 = np.radians(data['Euler3'].values)
    R    = euler_to_rotation_matrix(phi1, PHI, phi2)  # (N,3,3)

    pairs = grains._all_pairs
    if pairs is None or len(pairs) == 0:
        return np.full(n, np.nan)

    # Compute misorientation for all adjacent pairs (same grain only for KAM)
    R1 = R[pairs[:, 0]]
    R2 = R[pairs[:, 1]]
    M  = np.einsum('nij,nik->njk', R1, R2)
    tr = np.einsum('nii->n', M)
    ang_deg = np.degrees(np.arccos(np.clip((tr - 1) / 2, -1, 1)))

    same_grain = g[pairs[:, 0]] == g[pairs[:, 1]]
    use = same_grain & (ang_deg < max_angle_deg) & \
          (g[pairs[:, 0]] >= 0) & (g[pairs[:, 1]] >= 0)

    kam = np.full(n, np.nan)
    kam_sum = np.zeros(n)
    kam_cnt = np.zeros(n, dtype=np.int32)

    np.add.at(kam_sum, pairs[use, 0], ang_deg[use])
    np.add.at(kam_sum, pairs[use, 1], ang_deg[use])
    np.add.at(kam_cnt, pairs[use, 0], 1)
    np.add.at(kam_cnt, pairs[use, 1], 1)

    valid = (kam_cnt > 0) & (g >= 0)
    kam[valid] = kam_sum[valid] / kam_cnt[valid]
    return kam


def calc_grod(grains) -> np.ndarray:
    """
    Grain Reference Orientation Deviation (GROD) per pixel.

    Angle between each pixel orientation and its grain's mean orientation.

    Returns
    -------
    grod : (N,) float array (degrees), NaN for unindexed pixels
    """
    from .reconstruction import euler_to_rotation_matrix

    data = grains.data
    g    = grains.grain_ids
    n    = len(data)

    phi1_all = np.radians(data['Euler1'].values)
    PHI_all  = np.radians(data['Euler2'].values)
    phi2_all = np.radians(data['Euler3'].values)
    R_all    = euler_to_rotation_matrix(phi1_all, PHI_all, phi2_all)

    me = grains.mean_euler   # (n_grains, 3) degrees
    R_mean = euler_to_rotation_matrix(
        np.radians(me[:, 0]), np.radians(me[:, 1]), np.radians(me[:, 2])
    )

    id_to_idx = {int(gid): i for i, gid in enumerate(grains._valid_ids)}
    grod = np.full(n, np.nan)

    for i, gid in enumerate(grains._valid_ids):
        px = grains._pixel_idx[int(gid)]
        Rpx = R_all[px]            # (K,3,3)
        Rg  = R_mean[i]            # (3,3)
        M   = np.einsum('ij,kjl->kil', Rg.T, Rpx)
        tr  = np.einsum('kii->k', M)
        grod[px] = np.degrees(np.arccos(np.clip((tr - 1) / 2, -1, 1)))
    return grod


def calc_misorientation_distribution(grain_boundary,
                                      bins: int = 36) -> tuple[np.ndarray, np.ndarray]:
    """
    Length-weighted misorientation angle distribution.

    Parameters
    ----------
    grain_boundary : GrainBoundary object
    bins           : number of histogram bins in [0°, 180°]

    Returns
    -------
    centers : (bins,) bin centre angles (degrees)
    density : (bins,) length-weighted frequency
    """
    centers, counts = grain_boundary.angle_distribution(bins=bins)
    total = counts.sum()
    density = counts / total if total > 0 else counts
    return centers, density


def calc_texture_index(euler_angles_deg: np.ndarray,
                        sigma_deg: float = 10.0) -> float:
    """
    Estimate a simple texture index (J) from Euler angles using a Gaussian ODF kernel.

    J = ∫ f(g)² dg  ≈ kernel density auto-correlation.

    Parameters
    ----------
    euler_angles_deg : (N, 3) array of Euler angles (degrees)
    sigma_deg        : kernel half-width (degrees)

    Returns
    -------
    J : float ≥ 1 (1 = random texture)
    """
    if len(euler_angles_deg) == 0:
        return 1.0

    from .reconstruction import euler_to_rotation_matrix
    phi1 = np.radians(euler_angles_deg[:, 0])
    PHI  = np.radians(euler_angles_deg[:, 1])
    phi2 = np.radians(euler_angles_deg[:, 2])
    R = euler_to_rotation_matrix(phi1, PHI, phi2)   # (N,3,3)

    sigma = np.radians(sigma_deg)
    N = len(R)
    # Subsample for speed if large
    if N > 500:
        rng = np.random.default_rng(0)
        idx = rng.choice(N, 500, replace=False)
        R = R[idx]
        N = 500

    # Pairwise misorientation angles
    total = 0.0
    for i in range(N):
        M = np.einsum('ij,kjl->kil', R[i].T, R)
        tr = np.clip(np.einsum('kii->k', M), -1, 3)
        angles = np.arccos((tr - 1) / 2)
        total += np.sum(np.exp(-angles**2 / (2 * sigma**2)))

    J = total / N**2 / (sigma**3 * (2 * np.pi)**1.5) * 8 * np.pi**2
    return max(float(J), 1.0)


def subgrain_boundary_fraction(grains,
                                lagb_threshold: float = 15.0) -> np.ndarray:
    """
    Per-grain fraction of boundary length that is sub-grain (low-angle).

    Returns
    -------
    fractions : (n_grains,) float array in [0, 1]
    """
    gb = grains.boundary
    angles = gb.misorientation_angle
    lengths = gb.segment_length
    gp = gb.grain_id_pairs

    id_to_idx = {int(gid): i for i, gid in enumerate(grains._valid_ids)}
    total  = np.zeros(grains.n_grains)
    lagb   = np.zeros(grains.n_grains)

    for k, (ga, gb_id) in enumerate(gp):
        for gid in [ga, gb_id]:
            i = id_to_idx.get(int(gid))
            if i is not None:
                total[i] += lengths[k]
                if angles[k] < lagb_threshold:
                    lagb[i] += lengths[k]

    frac = np.where(total > 0, lagb / total, 0.0)
    return frac
