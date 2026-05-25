"""
Grain reconstruction from EBSD data using misorientation-threshold clustering.

Algorithm (equivalent to MTEX calcGrains with 'unitCell' option):
  1. Convert Euler angles to rotation matrices (vectorised)
  2. Build 4-connectivity pixel adjacency from XY grid
  3. Compute raw (no-symmetry) misorientation for every adjacent pair
  4. Pairs below the threshold (and same phase) remain connected
  5. scipy connected_components → grain labels
  6. Remove grains with fewer than min_pixels; label those pixels -1
"""

import numpy as np
import pandas as pd
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def euler_to_rotation_matrix(phi1: np.ndarray, PHI: np.ndarray,
                              phi2: np.ndarray) -> np.ndarray:
    """
    Vectorised Bunge ZXZ Euler angles (radians) → rotation matrices (N,3,3).

    R = Rz(φ1) @ Rx(Φ) @ Rz(φ2)
    """
    c1, s1 = np.cos(phi1), np.sin(phi1)
    cP, sP = np.cos(PHI),  np.sin(PHI)
    c2, s2 = np.cos(phi2), np.sin(phi2)

    n = len(phi1)
    R = np.empty((n, 3, 3), dtype=np.float64)
    R[:, 0, 0] =  c1*c2 - s1*cP*s2
    R[:, 0, 1] = -c1*s2 - s1*cP*c2
    R[:, 0, 2] =  s1*sP
    R[:, 1, 0] =  s1*c2 + c1*cP*s2
    R[:, 1, 1] = -s1*s2 + c1*cP*c2
    R[:, 1, 2] = -c1*sP
    R[:, 2, 0] =  sP*s2
    R[:, 2, 1] =  sP*c2
    R[:, 2, 2] =  cP
    return R


def build_pixel_adjacency(data: pd.DataFrame) -> np.ndarray:
    """
    Return (M, 2) array of pixel-index pairs connected in a 4-neighbourhood.

    Works with any regular square EBSD grid; missing pixels (holes) are handled
    by treating the grid as sparse.
    """
    x = data['X'].values.astype(np.float64)
    y = data['Y'].values.astype(np.float64)

    x_unique = np.sort(np.unique(x))
    y_unique = np.sort(np.unique(y))

    dx = np.min(np.diff(x_unique)) if len(x_unique) > 1 else 1.0
    dy = np.min(np.diff(y_unique)) if len(y_unique) > 1 else 1.0

    x_idx = np.round((x - x_unique[0]) / dx).astype(np.int32)
    y_idx = np.round((y - y_unique[0]) / dy).astype(np.int32)

    nx = int(x_idx.max()) + 1
    ny = int(y_idx.max()) + 1

    grid = np.full((ny, nx), -1, dtype=np.int64)
    grid[y_idx, x_idx] = np.arange(len(data), dtype=np.int64)

    pairs_list = []

    # Horizontal: grid[r, c] — grid[r, c+1]
    valid = (grid[:, :-1] >= 0) & (grid[:, 1:] >= 0)
    r, c = np.where(valid)
    pairs_list.append(np.column_stack([grid[r, c], grid[r, c + 1]]))

    # Vertical: grid[r, c] — grid[r+1, c]
    valid = (grid[:-1, :] >= 0) & (grid[1:, :] >= 0)
    r, c = np.where(valid)
    pairs_list.append(np.column_stack([grid[r, c], grid[r + 1, c]]))

    return np.vstack(pairs_list).astype(np.int64)


def raw_misorientation_deg(R: np.ndarray, pairs: np.ndarray) -> np.ndarray:
    """
    Unsymmetrised misorientation angles in degrees for each pixel pair.

    M = R1^T @ R2;  angle = arccos((tr(M)-1)/2)
    """
    R1 = R[pairs[:, 0]]
    R2 = R[pairs[:, 1]]
    M = np.einsum('nij,nik->njk', R1, R2)   # (N,3,3): R1^T @ R2
    trace = np.einsum('nii->n', M)
    return np.degrees(np.arccos(np.clip((trace - 1.0) / 2.0, -1.0, 1.0)))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def calc_grains(data: pd.DataFrame,
                threshold_deg: float = 10.0,
                min_pixels: int = 2) -> tuple[np.ndarray, np.ndarray]:
    """
    Segment an EBSD dataset into grains.

    Parameters
    ----------
    data          : DataFrame with columns X, Y, Phase, Euler1, Euler2, Euler3
    threshold_deg : misorientation threshold for grain boundary (degrees)
    min_pixels    : grains with fewer pixels are discarded (label → -1)

    Returns
    -------
    grain_ids : (N,) int array; -1 = unindexed / too small
    pairs     : (M, 2) adjacency pairs (useful for boundary detection)
    """
    n = len(data)
    phase = data['Phase'].values

    # Only indexed pixels (Phase > 0) participate in orientation clustering
    phi1 = np.radians(data['Euler1'].values)
    PHI  = np.radians(data['Euler2'].values)
    phi2 = np.radians(data['Euler3'].values)
    R = euler_to_rotation_matrix(phi1, PHI, phi2)

    pairs = build_pixel_adjacency(data)

    angles = raw_misorientation_deg(R, pairs)

    same_phase = phase[pairs[:, 0]] == phase[pairs[:, 1]]
    # Treat Phase == 0 as unindexed — always a boundary
    both_indexed = (phase[pairs[:, 0]] > 0) & (phase[pairs[:, 1]] > 0)
    keep = (angles < threshold_deg) & same_phase & both_indexed

    i = pairs[keep, 0]
    j = pairs[keep, 1]
    adj = csr_matrix((np.ones(len(i), dtype=np.float32), (i, j)),
                     shape=(n, n))
    adj = adj + adj.T

    _, labels = connected_components(adj, directed=False)

    # Relabel unindexed pixels
    labels = labels.astype(np.int32)
    labels[phase == 0] = -1

    # Remove small grains
    if min_pixels > 1:
        unique, counts = np.unique(labels[labels >= 0], return_counts=True)
        small = unique[counts < min_pixels]
        if len(small):
            mask = np.isin(labels, small)
            labels[mask] = -1

    # Re-number labels to be contiguous 0..K-1
    valid_ids = np.unique(labels[labels >= 0])
    remap = np.full(labels.max() + 2, -1, dtype=np.int32)
    remap[valid_ids] = np.arange(len(valid_ids), dtype=np.int32)
    grain_ids = np.where(labels >= 0, remap[labels], -1)

    return grain_ids, pairs
