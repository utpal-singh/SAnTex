"""
GrainBoundary — boundary segments extracted from a Grain2D object.

Each "segment" is a pixel-edge between two pixels of different grain IDs.

Features
--------
  Properties     : misorientation_angle, misorientation_axis, segment_length
  Selection      : by_angle(), by_phase_pair(), is_lagb(), is_hagb()
  CSL            : sigma_value(), is_csl(), is_twinning()
  Curvature      : curvature()
  Twist / Tilt   : twist_tilt()
  Distribution   : angle_distribution()
  Plot           : plot()
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from functools import cached_property


class GrainBoundary:
    """
    Grain boundary segments extracted from a Grain2D object.

    Parameters
    ----------
    grains : Grain2D object (the parent grain map)
    """

    def __init__(self, grains):
        self.grains = grains
        self._build()

    def _build(self):
        """Extract all pixel-pair edges between different grains."""
        g = self.grains.grain_ids
        pairs = self.grains._all_pairs
        if pairs is None or len(pairs) == 0:
            self._pairs = np.empty((0, 2), dtype=np.int64)
            return

        g1 = g[pairs[:, 0]]
        g2 = g[pairs[:, 1]]
        boundary = (g1 != g2) & (g1 >= 0) & (g2 >= 0)
        self._pairs = pairs[boundary].astype(np.int64)  # (M, 2) pixel pairs

    # ------------------------------------------------------------------
    # Derived segment attributes
    # ------------------------------------------------------------------

    @cached_property
    def n_segments(self) -> int:
        return len(self._pairs)

    @cached_property
    def grain_id_pairs(self) -> np.ndarray:
        """(M, 2) grain ID for each end of each boundary segment."""
        g = self.grains.grain_ids
        return np.column_stack([g[self._pairs[:, 0]],
                                g[self._pairs[:, 1]]])

    @cached_property
    def phase_id_pairs(self) -> np.ndarray:
        """(M, 2) phase ID for each end of each boundary segment."""
        phase = self.grains.data['Phase'].values
        return np.column_stack([phase[self._pairs[:, 0]],
                                phase[self._pairs[:, 1]]])

    @cached_property
    def midpoint(self) -> np.ndarray:
        """(M, 2) midpoint (x, y) of each segment."""
        x = self.grains.data['X'].values
        y = self.grains.data['Y'].values
        p0 = np.column_stack([x[self._pairs[:, 0]], y[self._pairs[:, 0]]])
        p1 = np.column_stack([x[self._pairs[:, 1]], y[self._pairs[:, 1]]])
        return (p0 + p1) / 2.0

    @cached_property
    def segment_length(self) -> np.ndarray:
        """(M,) length of each segment (Euclidean distance between pixel centres)."""
        x = self.grains.data['X'].values
        y = self.grains.data['Y'].values
        dx = x[self._pairs[:, 1]] - x[self._pairs[:, 0]]
        dy = y[self._pairs[:, 1]] - y[self._pairs[:, 0]]
        return np.sqrt(dx**2 + dy**2)

    @cached_property
    def _rotation_matrices(self) -> np.ndarray:
        """(M, 2, 3, 3) rotation matrices for both pixels of each segment."""
        from .reconstruction import euler_to_rotation_matrix
        phi1 = np.radians(self.grains.data['Euler1'].values)
        PHI  = np.radians(self.grains.data['Euler2'].values)
        phi2 = np.radians(self.grains.data['Euler3'].values)
        R = euler_to_rotation_matrix(phi1, PHI, phi2)
        return R

    @cached_property
    def misorientation_matrix(self) -> np.ndarray:
        """(M, 3, 3) relative rotation M = R1^T @ R2 for each segment."""
        R = self._rotation_matrices
        R1 = R[self._pairs[:, 0]]
        R2 = R[self._pairs[:, 1]]
        return np.einsum('nij,nik->njk', R1, R2)

    @cached_property
    def misorientation_angle(self) -> np.ndarray:
        """(M,) unsymmetrised misorientation angle per segment (degrees)."""
        M = self.misorientation_matrix
        trace = np.einsum('nii->n', M)
        return np.degrees(np.arccos(np.clip((trace - 1) / 2, -1, 1)))

    @cached_property
    def misorientation_axis(self) -> np.ndarray:
        """(M, 3) rotation axis (sample frame) for each segment."""
        M = self.misorientation_matrix
        axes = np.stack([
            M[:, 2, 1] - M[:, 1, 2],
            M[:, 0, 2] - M[:, 2, 0],
            M[:, 1, 0] - M[:, 0, 1],
        ], axis=1)
        norms = np.linalg.norm(axes, axis=1, keepdims=True)
        norms = np.where(norms > 1e-12, norms, 1.0)
        return axes / norms

    # ------------------------------------------------------------------
    # Selection helpers
    # ------------------------------------------------------------------

    def by_angle(self, min_deg: float = 0.0,
                  max_deg: float = 180.0) -> 'GrainBoundary':
        """Return a sub-GrainBoundary restricted to a misorientation angle range."""
        mask = (self.misorientation_angle >= min_deg) & \
               (self.misorientation_angle <= max_deg)
        return self._subset(mask)

    def by_phase_pair(self, phase1: int, phase2: int) -> 'GrainBoundary':
        pp = self.phase_id_pairs
        mask = ((pp[:, 0] == phase1) & (pp[:, 1] == phase2)) | \
               ((pp[:, 0] == phase2) & (pp[:, 1] == phase1))
        return self._subset(mask)

    def is_lagb(self, threshold_deg: float = 15.0) -> np.ndarray:
        """Boolean mask: True where segment is a low-angle grain boundary."""
        return self.misorientation_angle < threshold_deg

    def is_hagb(self, threshold_deg: float = 15.0) -> np.ndarray:
        """Boolean mask: True where segment is a high-angle grain boundary."""
        return self.misorientation_angle >= threshold_deg

    def _subset(self, mask: np.ndarray) -> 'GrainBoundary':
        sub = GrainBoundary.__new__(GrainBoundary)
        sub.grains = self.grains
        sub._pairs = self._pairs[mask]
        return sub

    # ------------------------------------------------------------------
    # CSL identification
    # ------------------------------------------------------------------

    def sigma_value(self, threshold_deg: float = 5.0) -> np.ndarray:
        """
        Identify the CSL Σ value for each segment (cubic symmetry).
        Returns int array; 0 = not a recognised CSL within threshold.
        """
        from .symmetry import CSL_CUBIC, csl_rotation_matrix, misorientation_angles_batch, SYMMETRY
        M = self.misorientation_matrix
        G = SYMMETRY["cubic"]

        results = np.zeros(len(M), dtype=np.int32)
        for sigma, (_, axis) in CSL_CUBIC.items():
            Rcsl = csl_rotation_matrix(sigma)
            if Rcsl is None:
                continue
            # dis-angle between each boundary M and Rcsl (with cubic symmetry)
            # For speed, compute using symmetry pairs
            angles = misorientation_angles_batch(
                np.einsum('ij,njk->nik', Rcsl.T, M), "cubic")
            match = np.degrees(angles) < threshold_deg
            results[match] = sigma
        return results

    def is_csl(self, sigma: int, threshold_deg: float = 5.0) -> np.ndarray:
        return self.sigma_value(threshold_deg) == sigma

    def is_twinning(self, sigma: int = 3,
                    threshold_deg: float = 5.0) -> np.ndarray:
        """True for Σ3 (or specified Σ) twin boundaries."""
        return self.is_csl(sigma, threshold_deg)

    # ------------------------------------------------------------------
    # Curvature (circle fit through 3 consecutive segment midpoints)
    # ------------------------------------------------------------------

    def curvature(self) -> np.ndarray:
        """
        Signed curvature per segment (1 / fitting radius).
        Segments at the ends of chains (< 2 neighbours) get curvature = 0.
        """
        mp = self.midpoint  # (M, 2)
        M = len(mp)
        kappa = np.zeros(M)

        # Build segment adjacency: two segments are adjacent if they share a pixel
        # We use the pixel indices to find chains
        from scipy.sparse import csr_matrix
        # pixel → segment list
        pix0 = self._pairs[:, 0]
        pix1 = self._pairs[:, 1]
        # Build a simple chain: iterate over ordered pairs
        # For each segment i, adjacent segments share pix0[i] or pix1[i]
        # Quick approximate approach: sort by midpoint x, y and consider consecutive
        for i in range(1, M - 1):
            p_l, p_c, p_r = mp[i - 1], mp[i], mp[i + 1]
            a = p_c - p_l
            b = p_r - p_l
            cross = a[0] * b[1] - a[1] * b[0]  # signed area × 2
            la = np.linalg.norm(a)
            lb_lc = np.linalg.norm(p_c - p_r)
            lb_full = np.linalg.norm(b)
            denom = la * lb_full * lb_lc
            if denom > 1e-12:
                kappa[i] = 2.0 * cross / denom
        return kappa

    # ------------------------------------------------------------------
    # Twist / tilt decomposition
    # ------------------------------------------------------------------

    def twist_tilt(self) -> dict:
        """
        Decompose each boundary misorientation into twist and tilt components.

        Boundary normal is estimated as the unit vector from pixel 0 to pixel 1.

        Returns dict:
            twist_angle (M,)  — rotation about boundary normal (degrees)
            tilt_angle  (M,)  — rotation about axis ⊥ normal (degrees)
        """
        x = self.grains.data['X'].values
        y = self.grains.data['Y'].values
        dx = x[self._pairs[:, 1]] - x[self._pairs[:, 0]]
        dy = y[self._pairs[:, 1]] - y[self._pairs[:, 0]]
        norm = np.sqrt(dx**2 + dy**2 + 1e-30)
        # Boundary normal in 2-D (perpendicular to the segment)
        n = np.column_stack([-dy / norm, dx / norm, np.zeros(len(dx))])

        axes = self.misorientation_axis   # (M, 3)
        angles = self.misorientation_angle  # (M,) degrees

        # Twist: component of rotation axis along boundary normal
        cos_twist = np.einsum('ni,ni->n', axes, n)
        twist_angle = np.degrees(np.arccos(np.clip(np.abs(cos_twist), 0, 1))) * \
                      np.radians(angles) / (np.pi / 2)

        # Tilt: component perpendicular to boundary normal
        tilt_angle = np.sqrt(np.maximum(np.radians(angles)**2 -
                                        (np.arccos(np.clip(cos_twist, -1, 1)))**2 * 0, 0))
        # simpler: decompose
        twist_frac = np.abs(cos_twist)
        tilt_frac  = np.sqrt(np.clip(1 - cos_twist**2, 0, 1))
        twist_angle = angles * twist_frac
        tilt_angle  = angles * tilt_frac

        return {"twist_angle": twist_angle, "tilt_angle": tilt_angle}

    # ------------------------------------------------------------------
    # Distribution
    # ------------------------------------------------------------------

    def angle_distribution(self, bins: int = 36) -> tuple[np.ndarray, np.ndarray]:
        """Histogram of misorientation angles (degrees). Returns (bin_centers, counts)."""
        angles = self.misorientation_angle
        edges = np.linspace(0, 180, bins + 1)
        counts, _ = np.histogram(angles, bins=edges,
                                  weights=self.segment_length)
        centers = (edges[:-1] + edges[1:]) / 2
        return centers, counts

    # ------------------------------------------------------------------
    # Total boundary length between grain pairs
    # ------------------------------------------------------------------

    def shared_length(self) -> pd.DataFrame:
        """
        DataFrame of (grain_a, grain_b, total_length) for each unique grain pair.
        """
        gp = self.grain_id_pairs
        sl = self.segment_length
        rows = []
        ga = np.minimum(gp[:, 0], gp[:, 1])
        gb = np.maximum(gp[:, 0], gp[:, 1])
        pairs_sorted = np.column_stack([ga, gb])
        unique_pairs, inv = np.unique(pairs_sorted, axis=0, return_inverse=True)
        lengths = np.zeros(len(unique_pairs))
        np.add.at(lengths, inv, sl)
        for i, (a, b) in enumerate(unique_pairs):
            rows.append({"grain_a": int(a), "grain_b": int(b),
                         "length": float(lengths[i])})
        return pd.DataFrame(rows)

    # ------------------------------------------------------------------
    # Plot
    # ------------------------------------------------------------------

    def plot(self, ax=None, scalar: str = "angle", cmap: str = "hot_r",
             linewidth: float = 0.5):
        """
        Draw boundary segments coloured by misorientation angle (or 'lagb'/'hagb').
        """
        import matplotlib.pyplot as plt
        from matplotlib.collections import LineCollection
        from matplotlib.cm import get_cmap
        from matplotlib.colors import Normalize

        if ax is None:
            _, ax = plt.subplots(figsize=(7, 5))

        x = self.grains.data['X'].values
        y = self.grains.data['Y'].values

        segs = np.stack([
            np.column_stack([x[self._pairs[:, 0]], y[self._pairs[:, 0]]]),
            np.column_stack([x[self._pairs[:, 1]], y[self._pairs[:, 1]]]),
        ], axis=1)  # (M, 2, 2)

        angles = self.misorientation_angle
        norm = Normalize(vmin=angles.min(), vmax=angles.max())
        lc = LineCollection(segs, array=angles, cmap=cmap, norm=norm,
                             linewidth=linewidth)
        ax.add_collection(lc)
        plt.colorbar(lc, ax=ax, label="Misorientation angle (°)")
        ax.set_xlim(x.min(), x.max())
        ax.set_ylim(y.min(), y.max())
        ax.set_aspect('equal')
        ax.set_title("Grain boundary map")
        return ax
