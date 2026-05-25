"""
calcGrainBoundary.py — grain reconstruction helpers for santex.ebsd.

The original implementation had six fundamental errors (documented below).
This file provides corrected drop-in replacements that match the MTEX algorithm.

Original bugs
-------------
1. Wrong algorithm
   The original "compare against grain representative" sequential loop is a
   naive orientation-clustering algorithm that ignores spatial adjacency.
   Two pixels at opposite corners of the map can end up in the same grain if
   they happen to have similar orientations.
   FIX → build 4-connectivity pixel adjacency from XY grid, threshold each
   adjacent pixel pair, run connected-components (scipy).  Matches MTEX
   doSegmentation() / connectedComponents().

2. Wrong quaternion formula (Euler → quaternion)
   The original formula produces a rotation in an unidentified convention that
   is NOT Bunge ZXZ. Tested: euler_to_quaternion(phi=π/2, 0, 0) returns a 90°
   rotation about the X-axis, not Z.
   FIX → use the standard Bunge ZXZ formula:
     q0 = cos(Φ/2)·cos((φ1+φ2)/2)  = c1·c2·c3 − s1·c2·s3
     q1 = sin(Φ/2)·cos((φ1−φ2)/2)  = c1·s2·c3 + s1·s2·s3
     q2 = sin(Φ/2)·sin((φ1−φ2)/2)  = s1·s2·c3 − c1·s2·s3
     q3 = cos(Φ/2)·sin((φ1+φ2)/2)  = s1·c2·c3 + c1·c2·s3
   where c1=cos(φ1/2), s1=sin(φ1/2), etc.

3. No crystal symmetry reduction
   MTEX gbc_angle.m: d = max|dot(mori, cs.properGroup.rot)|; criterion based on
   cos(threshold/2). Without symmetry reduction the raw quaternion angle can
   vastly overestimate the true disorientation.
   FIX → apply symmetry operators and take the minimum-angle (disorientation).

4. Wrong grain mean orientation update
   np.mean([grain_euler, pixel_euler]) is an arithmetic mean of Euler angles.
   This is geometrically meaningless (singularities, 2π wrap-around).
   FIX → quaternion averaging via rotation-matrix mean + SVD re-orthogonalisation.

5. Race condition in parallel implementation
   Parallel(prefer="threads") on shared mutable lists causes non-deterministic
   results (unsynchronised reads and writes).
   FIX → no parallel mutable state; the new algorithm is fully vectorised numpy.

6. Non-contiguous grains
   Without spatial adjacency, distant pixels with similar orientations are
   merged into one "grain", violating the definition of a grain.
   FIX → only pixels sharing an edge in the XY grid can belong to the same grain.
"""

from __future__ import annotations
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Corrected Bunge ZXZ Euler → quaternion
# ---------------------------------------------------------------------------

def euler_to_quaternion(phi1: float, PHI: float, phi2: float) -> np.ndarray:
    """
    Convert Bunge ZXZ Euler angles (radians) to a unit quaternion [w, x, y, z].

    R = Rz(φ1) @ Rx(Φ) @ Rz(φ2)

    Parameters
    ----------
    phi1, PHI, phi2 : Euler angles in radians (floats or arrays)

    Returns
    -------
    q : shape (..., 4) unit quaternion [w, x, y, z]
    """
    phi1 = np.asarray(phi1, dtype=np.float64)
    PHI  = np.asarray(PHI,  dtype=np.float64)
    phi2 = np.asarray(phi2, dtype=np.float64)

    c1, s1 = np.cos(phi1 / 2), np.sin(phi1 / 2)
    c2, s2 = np.cos(PHI  / 2), np.sin(PHI  / 2)
    c3, s3 = np.cos(phi2 / 2), np.sin(phi2 / 2)

    w = c1 * c2 * c3 - s1 * c2 * s3
    x = c1 * s2 * c3 + s1 * s2 * s3
    y = s1 * s2 * c3 - c1 * s2 * s3
    z = s1 * c2 * c3 + c1 * c2 * s3

    return np.stack([w, x, y, z], axis=-1)


def euler_to_quaternion_deg(phi1_deg, PHI_deg, phi2_deg) -> np.ndarray:
    """Convenience wrapper — inputs in degrees."""
    return euler_to_quaternion(
        np.radians(phi1_deg), np.radians(PHI_deg), np.radians(phi2_deg))


# ---------------------------------------------------------------------------
# Symmetry-reduced misorientation (disorientation)
# ---------------------------------------------------------------------------

def misorientation_angle(euler1: np.ndarray, euler2: np.ndarray,
                          symmetry: str = "triclinic") -> float:
    """
    Compute the disorientation angle (degrees) between two orientations.

    Parameters
    ----------
    euler1, euler2 : array-like of shape (3,), Euler angles in degrees
    symmetry       : crystal symmetry ('cubic', 'hexagonal', 'orthorhombic', etc.)

    Returns
    -------
    angle : disorientation angle in degrees (minimum over crystal symmetry)
    """
    q1 = euler_to_quaternion(*np.radians(euler1))
    q2 = euler_to_quaternion(*np.radians(euler2))
    return _disorientation_deg(q1, q2, symmetry)


def _disorientation_deg(q1: np.ndarray, q2: np.ndarray,
                         symmetry: str = "triclinic") -> float:
    """
    Disorientation angle between two unit quaternions [w,x,y,z].

    Matches MTEX gbc_angle: d = max|dot(mori, sym_ops)|; angle = 2*arccos(d)
    """
    from ..grains.symmetry import SYMMETRY
    G = SYMMETRY.get(symmetry.lower(), SYMMETRY["triclinic"])

    # Relative quaternion: q_rel = inv(q1) * q2
    # inv(q) = conjugate for unit quaternion = [w, -x, -y, -z]
    q1_inv = q1 * np.array([1, -1, -1, -1])
    q_rel = _quat_mul(q1_inv, q2)

    # Convert symmetry rotation matrices to quaternions
    max_dot = 0.0
    for Rg in G:
        q_sym = _rot_to_quat(Rg)
        # dot product (equivalent to cos(angle/2) of relative rotation)
        dot_val = abs(np.dot(q_rel, q_sym))
        if dot_val > max_dot:
            max_dot = dot_val

    max_dot = min(max_dot, 1.0)
    return float(np.degrees(2.0 * np.arccos(max_dot)))


def _quat_mul(q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
    """Quaternion multiplication q1 * q2, both [w, x, y, z]."""
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2,
    ])


def _rot_to_quat(R: np.ndarray) -> np.ndarray:
    """Convert a 3×3 rotation matrix to a unit quaternion [w, x, y, z]."""
    trace = R[0, 0] + R[1, 1] + R[2, 2]
    if trace > 0:
        s = 0.5 / np.sqrt(trace + 1.0)
        w = 0.25 / s
        x = (R[2, 1] - R[1, 2]) * s
        y = (R[0, 2] - R[2, 0]) * s
        z = (R[1, 0] - R[0, 1]) * s
    elif R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
        s = 2.0 * np.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2])
        w = (R[2, 1] - R[1, 2]) / s
        x = 0.25 * s
        y = (R[0, 1] + R[1, 0]) / s
        z = (R[0, 2] + R[2, 0]) / s
    elif R[1, 1] > R[2, 2]:
        s = 2.0 * np.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2])
        w = (R[0, 2] - R[2, 0]) / s
        x = (R[0, 1] + R[1, 0]) / s
        y = 0.25 * s
        z = (R[1, 2] + R[2, 1]) / s
    else:
        s = 2.0 * np.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1])
        w = (R[1, 0] - R[0, 1]) / s
        x = (R[0, 2] + R[2, 0]) / s
        y = (R[1, 2] + R[2, 1]) / s
        z = 0.25 * s
    q = np.array([w, x, y, z])
    return q / np.linalg.norm(q)


# ---------------------------------------------------------------------------
# Corrected grain assignment — spatial adjacency + connected components
# ---------------------------------------------------------------------------

def assign_to_grains(data: pd.DataFrame,
                     threshold: float = 10.0,
                     min_pixels: int = 2) -> dict[int, int]:
    """
    Assign grain IDs to each pixel via spatial adjacency + misorientation threshold.

    This replaces the buggy `assign_to_grains_parallel`.  It is fully
    vectorised (no Python loops over pixels) and produces spatially contiguous
    grains, matching the MTEX calcGrains algorithm.

    Parameters
    ----------
    data      : DataFrame with columns X, Y, Phase, Euler1, Euler2, Euler3
    threshold : grain boundary misorientation angle threshold (degrees)
    min_pixels: discard grains smaller than this (label → -1)

    Returns
    -------
    grain_indices : dict mapping DataFrame row index → grain ID (int)
                   grain ID = -1 for unindexed / too-small pixels
    """
    from ..grains.reconstruction import calc_grains
    grain_ids, _ = calc_grains(data.reset_index(drop=True), threshold, min_pixels)
    return dict(enumerate(grain_ids.tolist()))


# Backwards-compatible alias (original name)
def assign_to_grains_parallel(df: pd.DataFrame,
                               threshold: float = 10.0) -> dict[int, int]:
    """
    Drop-in replacement for the original assign_to_grains_parallel.

    The 'parallel' in the name is a historical artefact; the new implementation
    is vectorised numpy/scipy and does not require multi-threading.

    See assign_to_grains() for full documentation.
    """
    return assign_to_grains(df, threshold=threshold, min_pixels=2)
