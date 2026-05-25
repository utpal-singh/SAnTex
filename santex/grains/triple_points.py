"""
TriplePoints — triple junctions in a Grain2D map.

A triple point is a pixel that is adjacent to pixels from (at least) 3 distinct grains.
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from functools import cached_property


class TriplePoints:
    """
    Triple junctions extracted from a Grain2D object.

    Detection: scan every pixel; if its 4-neighbours contain pixels from ≥ 3
    distinct grains (counting the pixel itself), it is a triple point.
    """

    def __init__(self, grains):
        self.grains = grains
        self._detect()

    def _detect(self):
        from .reconstruction import build_pixel_adjacency

        g = self.grains.grain_ids
        pairs = self.grains._all_pairs

        if pairs is None or len(pairs) == 0:
            self._tp_idx = np.array([], dtype=np.int64)
            self._tp_grain_ids = np.empty((0, 4), dtype=np.int32)
            return

        # Build neighbour lists efficiently using adjacency pairs
        n = len(g)
        from collections import defaultdict
        neighbours: dict[int, set] = defaultdict(set)
        for i, j in pairs:
            neighbours[int(i)].add(int(j))
            neighbours[int(j)].add(int(i))

        tp_idx = []
        tp_grains = []  # list of sets of grain IDs meeting at each tp

        for px, nbs in neighbours.items():
            if g[px] < 0:
                continue
            grain_set = {int(g[px])}
            for nb in nbs:
                if g[nb] >= 0:
                    grain_set.add(int(g[nb]))
            if len(grain_set) >= 3:
                tp_idx.append(px)
                # Store the 3 most common grain IDs (up to 4)
                sorted_gs = sorted(grain_set)[:4]
                while len(sorted_gs) < 4:
                    sorted_gs.append(-1)
                tp_grains.append(sorted_gs[:4])

        self._tp_idx = np.array(tp_idx, dtype=np.int64)
        if tp_grains:
            self._tp_grain_ids = np.array(tp_grains, dtype=np.int32)
        else:
            self._tp_grain_ids = np.empty((0, 4), dtype=np.int32)

    # ------------------------------------------------------------------

    @property
    def n(self) -> int:
        return len(self._tp_idx)

    @cached_property
    def coordinates(self) -> np.ndarray:
        """(N, 2) x, y coordinates of each triple point."""
        x = self.grains.data['X'].values
        y = self.grains.data['Y'].values
        return np.column_stack([x[self._tp_idx], y[self._tp_idx]])

    @property
    def grain_ids(self) -> np.ndarray:
        """(N, 4) grain IDs meeting at each triple point (-1 = padding)."""
        return self._tp_grain_ids

    @cached_property
    def misorientation_angles(self) -> np.ndarray:
        """
        (N, 3) pairwise misorientation angles (degrees) between the first 3
        grains meeting at each triple point.
        """
        from .reconstruction import euler_to_rotation_matrix
        phi1 = np.radians(self.grains.data['Euler1'].values)
        PHI  = np.radians(self.grains.data['Euler2'].values)
        phi2 = np.radians(self.grains.data['Euler3'].values)
        R = euler_to_rotation_matrix(phi1, PHI, phi2)

        id_to_idx = {int(gid): i for i, gid in enumerate(self.grains._valid_ids)}
        # mean rotation per grain
        from .reconstruction import euler_to_rotation_matrix as e2r
        me = self.grains.mean_euler  # (n_grains, 3)
        R_mean = e2r(np.radians(me[:, 0]), np.radians(me[:, 1]), np.radians(me[:, 2]))

        out = np.zeros((self.n, 3), dtype=np.float64)
        for k in range(self.n):
            gs = self._tp_grain_ids[k]
            valid = [g for g in gs if g >= 0 and g in id_to_idx]
            if len(valid) < 2:
                continue
            pairs_tp = [(valid[0], valid[1]),
                        (valid[0], valid[2] if len(valid) > 2 else valid[1]),
                        (valid[1], valid[2] if len(valid) > 2 else valid[1])]
            for p_idx, (ga, gb) in enumerate(pairs_tp[:3]):
                ia, ib = id_to_idx.get(ga), id_to_idx.get(gb)
                if ia is None or ib is None:
                    continue
                M = R_mean[ia].T @ R_mean[ib]
                trace = np.clip(np.trace(M), -1, 3)
                out[k, p_idx] = np.degrees(np.arccos((trace - 1) / 2))
        return out

    def to_dataframe(self) -> pd.DataFrame:
        coords = self.coordinates
        mori = self.misorientation_angles
        df = pd.DataFrame({
            "x": coords[:, 0],
            "y": coords[:, 1],
            "grain_a": self._tp_grain_ids[:, 0],
            "grain_b": self._tp_grain_ids[:, 1],
            "grain_c": self._tp_grain_ids[:, 2],
            "mori_ab_deg": mori[:, 0],
            "mori_ac_deg": mori[:, 1],
            "mori_bc_deg": mori[:, 2],
        })
        return df

    def plot(self, ax=None, color: str = "red", marker: str = "x",
             markersize: float = 4):
        import matplotlib.pyplot as plt
        if ax is None:
            _, ax = plt.subplots(figsize=(7, 5))
        if self.n > 0:
            coords = self.coordinates
            ax.plot(coords[:, 0], coords[:, 1], marker, color=color,
                    markersize=markersize, linewidth=0, label="Triple points")
            ax.legend(fontsize=8)
        ax.set_aspect('equal')
        return ax
