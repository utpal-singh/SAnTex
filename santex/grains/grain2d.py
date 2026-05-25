"""
Grain2D — main class representing a segmented EBSD dataset.

Provides:
  Shape:        area, perimeter, equivalent_radius, aspect_ratio, diameter
  Ellipse fit:  fit_ellipse()
  Convex hull:  convex_hull_props()
  Caliper:      feret_diameter()
  Orientation:  mean_euler, gos
  Topology:     neighbors(), merge()
  Selection:    select_by_phase(), select_by_size(), select_by_orientation()
  Plot:         plot_phase_map(), plot_grain_size_map()
  Summary:      summary_df()
"""

from __future__ import annotations
import numpy as np
import pandas as pd
from functools import cached_property
from scipy.spatial import ConvexHull
from scipy.spatial.distance import pdist

from .reconstruction import euler_to_rotation_matrix, calc_grains, build_pixel_adjacency, raw_misorientation_deg


class Grain2D:
    """
    Container for grain-segmented EBSD data.

    Parameters
    ----------
    data       : EBSD DataFrame (X, Y, Phase, Euler1, Euler2, Euler3)
    grain_ids  : (N,) int array from calc_grains(); -1 = unindexed
    all_pairs  : (M, 2) full pixel adjacency (output of calc_grains)
    """

    def __init__(self, data: pd.DataFrame, grain_ids: np.ndarray,
                 all_pairs: np.ndarray | None = None):
        self.data = data.reset_index(drop=True)
        self.grain_ids = grain_ids.astype(np.int32)
        self._all_pairs = all_pairs

        self._valid_ids = np.unique(grain_ids[grain_ids >= 0])
        self.n_grains = len(self._valid_ids)

        # pixel-index sets per grain (dict)
        self._pixel_idx: dict[int, np.ndarray] = {
            int(gid): np.where(grain_ids == gid)[0]
            for gid in self._valid_ids
        }

    # ------------------------------------------------------------------
    # Grid geometry
    # ------------------------------------------------------------------

    @cached_property
    def _step(self) -> tuple[float, float]:
        x_u = np.sort(np.unique(self.data['X'].values))
        y_u = np.sort(np.unique(self.data['Y'].values))
        dx = float(np.min(np.diff(x_u))) if len(x_u) > 1 else 1.0
        dy = float(np.min(np.diff(y_u))) if len(y_u) > 1 else 1.0
        return dx, dy

    @property
    def step_x(self) -> float:
        return self._step[0]

    @property
    def step_y(self) -> float:
        return self._step[1]

    # ------------------------------------------------------------------
    # Basic scalar properties (per grain, same order as _valid_ids)
    # ------------------------------------------------------------------

    @cached_property
    def grain_id_list(self) -> np.ndarray:
        return self._valid_ids.copy()

    @cached_property
    def phase_id(self) -> np.ndarray:
        """Dominant phase per grain (mode of pixel phases)."""
        phase = self.data['Phase'].values
        result = []
        for gid in self._valid_ids:
            px = self._pixel_idx[int(gid)]
            phases, counts = np.unique(phase[px], return_counts=True)
            result.append(int(phases[np.argmax(counts)]))
        return np.array(result, dtype=np.int32)

    @cached_property
    def n_pixels(self) -> np.ndarray:
        return np.array([len(self._pixel_idx[int(g)]) for g in self._valid_ids],
                        dtype=np.int64)

    @cached_property
    def area(self) -> np.ndarray:
        """Grain area in map-coordinate units²."""
        return self.n_pixels * self.step_x * self.step_y

    @cached_property
    def centroid(self) -> np.ndarray:
        """(n_grains, 2) array of (x, y) centroids."""
        x = self.data['X'].values
        y = self.data['Y'].values
        out = []
        for gid in self._valid_ids:
            px = self._pixel_idx[int(gid)]
            out.append([x[px].mean(), y[px].mean()])
        return np.array(out, dtype=np.float64)

    @cached_property
    def equivalent_radius(self) -> np.ndarray:
        """Radius of circle with same area as grain."""
        return np.sqrt(self.area / np.pi)

    @cached_property
    def equivalent_perimeter(self) -> np.ndarray:
        return 2.0 * np.pi * self.equivalent_radius

    @cached_property
    def perimeter(self) -> np.ndarray:
        """
        Estimated perimeter: count boundary pixel edges × step size.
        A boundary edge is an edge between pixels of different grain IDs.
        """
        if self._all_pairs is None:
            return np.zeros(self.n_grains)
        g = self.grain_ids
        pairs = self._all_pairs
        # edges where the two pixels belong to different grains (or one is -1)
        g1 = g[pairs[:, 0]]
        g2 = g[pairs[:, 1]]
        boundary_mask = g1 != g2
        b_pairs = pairs[boundary_mask]

        # count boundary edges per grain
        step = (self.step_x + self.step_y) / 2.0
        perim = np.zeros(self.n_grains, dtype=np.float64)
        for idx, gid in enumerate(self._valid_ids):
            n_edges = np.sum((g[b_pairs[:, 0]] == gid) |
                             (g[b_pairs[:, 1]] == gid))
            perim[idx] = n_edges * step
        return perim

    @cached_property
    def diameter(self) -> np.ndarray:
        """
        Maximum Feret diameter (longest caliper) per grain.
        Computed on the convex hull vertices for efficiency.
        """
        x = self.data['X'].values
        y = self.data['Y'].values
        diameters = []
        for gid in self._valid_ids:
            px = self._pixel_idx[int(gid)]
            pts = np.column_stack([x[px], y[px]])
            if len(pts) < 2:
                diameters.append(0.0)
                continue
            if len(pts) <= 3:
                diameters.append(float(pdist(pts).max()))
                continue
            try:
                hull = ConvexHull(pts)
                hull_pts = pts[hull.vertices]
                diameters.append(float(pdist(hull_pts).max()))
            except Exception:
                diameters.append(float(pdist(pts).max()))
        return np.array(diameters, dtype=np.float64)

    # ------------------------------------------------------------------
    # Ellipse fitting (PCA on boundary pixels)
    # ------------------------------------------------------------------

    def fit_ellipse(self) -> dict:
        """
        Fit an ellipse to each grain using PCA of pixel positions.

        Returns dict with arrays (all length n_grains):
            a, b    — semi-major / semi-minor axes (map units)
            angle   — orientation of major axis (degrees, CCW from +x)
            cx, cy  — centroid
        """
        x = self.data['X'].values
        y = self.data['Y'].values
        a_list, b_list, angle_list = [], [], []
        for gid in self._valid_ids:
            px = self._pixel_idx[int(gid)]
            pts = np.column_stack([x[px], y[px]])
            cx, cy = pts[:, 0].mean(), pts[:, 1].mean()
            v = pts - np.array([cx, cy])
            cov = (v.T @ v) / max(len(v) - 1, 1)
            eigvals, eigvecs = np.linalg.eigh(cov)
            eigvals = np.clip(eigvals, 0, None)
            # scale axes so ellipse area ≈ grain area
            scale = np.sqrt(self.step_x * self.step_y * np.pi)
            a = np.sqrt(eigvals[1]) * 2.0
            b = np.sqrt(eigvals[0]) * 2.0
            angle = np.degrees(np.arctan2(eigvecs[1, 1], eigvecs[0, 1]))
            a_list.append(a)
            b_list.append(b)
            angle_list.append(angle)
        return {
            "a":     np.array(a_list),
            "b":     np.array(b_list),
            "angle": np.array(angle_list),
            "cx":    self.centroid[:, 0],
            "cy":    self.centroid[:, 1],
        }

    @cached_property
    def aspect_ratio(self) -> np.ndarray:
        """Ratio a/b from ellipse fit (≥ 1)."""
        ell = self.fit_ellipse()
        b = np.where(ell["b"] > 0, ell["b"], 1e-12)
        return np.maximum(ell["a"] / b, 1.0)

    # ------------------------------------------------------------------
    # Convex hull properties
    # ------------------------------------------------------------------

    def convex_hull_props(self) -> dict:
        """
        Returns dict with arrays hull_area, hull_perimeter, solidity (n_grains).
        Solidity = grain_area / hull_area.
        """
        x = self.data['X'].values
        y = self.data['Y'].values
        hull_area, hull_perim, solidity = [], [], []
        for idx, gid in enumerate(self._valid_ids):
            px = self._pixel_idx[int(gid)]
            pts = np.column_stack([x[px], y[px]])
            if len(pts) < 4:
                hull_area.append(self.area[idx])
                hull_perim.append(self.perimeter[idx])
                solidity.append(1.0)
                continue
            try:
                hull = ConvexHull(pts)
                hv = pts[hull.vertices]
                # Shoelace formula for hull area
                ha = 0.5 * abs(np.dot(hv[:, 0], np.roll(hv[:, 1], -1)) -
                               np.dot(hv[:, 1], np.roll(hv[:, 0], -1)))
                # Hull perimeter
                hp = float(np.sum(np.linalg.norm(np.diff(hv, axis=0, append=hv[:1]), axis=1)))
                hull_area.append(ha)
                hull_perim.append(hp)
                solidity.append(self.area[idx] / ha if ha > 0 else 1.0)
            except Exception:
                hull_area.append(self.area[idx])
                hull_perim.append(self.perimeter[idx])
                solidity.append(1.0)
        return {
            "hull_area":      np.array(hull_area),
            "hull_perimeter": np.array(hull_perim),
            "solidity":       np.array(solidity),
        }

    # ------------------------------------------------------------------
    # Caliper / Feret diameter (rotating calipers, 72 directions)
    # ------------------------------------------------------------------

    def feret_diameter(self, n_angles: int = 72) -> dict:
        """
        Max and min Feret (caliper) diameters per grain.

        Returns dict: max_feret, min_feret, max_feret_angle, min_feret_angle
        """
        x = self.data['X'].values
        y = self.data['Y'].values
        angles = np.linspace(0, np.pi, n_angles, endpoint=False)
        dirs = np.column_stack([np.cos(angles), np.sin(angles)])  # (n_angles, 2)

        max_f, min_f, max_a, min_a = [], [], [], []
        for gid in self._valid_ids:
            px = self._pixel_idx[int(gid)]
            pts = np.column_stack([x[px], y[px]])
            if len(pts) < 2:
                max_f.append(0.0); min_f.append(0.0)
                max_a.append(0.0); min_a.append(0.0)
                continue
            try:
                hull = ConvexHull(pts)
                hpts = pts[hull.vertices]
            except Exception:
                hpts = pts
            proj = hpts @ dirs.T      # (n_hull, n_angles)
            widths = proj.max(axis=0) - proj.min(axis=0)  # (n_angles,)
            imax = int(np.argmax(widths))
            imin = int(np.argmin(widths))
            max_f.append(float(widths[imax]))
            min_f.append(float(widths[imin]))
            max_a.append(float(np.degrees(angles[imax])))
            min_a.append(float(np.degrees(angles[imin])))
        return {
            "max_feret":       np.array(max_f),
            "min_feret":       np.array(min_f),
            "max_feret_angle": np.array(max_a),
            "min_feret_angle": np.array(min_a),
        }

    # ------------------------------------------------------------------
    # Orientation parameters
    # ------------------------------------------------------------------

    @cached_property
    def mean_euler(self) -> np.ndarray:
        """
        Mean Euler angles per grain (degrees), shape (n_grains, 3).
        Uses rotation-matrix averaging (Frobenius-norm closest orthogonal matrix).
        """
        phi1 = np.radians(self.data['Euler1'].values)
        PHI  = np.radians(self.data['Euler2'].values)
        phi2 = np.radians(self.data['Euler3'].values)
        R_all = euler_to_rotation_matrix(phi1, PHI, phi2)

        out = np.zeros((self.n_grains, 3), dtype=np.float64)
        for idx, gid in enumerate(self._valid_ids):
            px = self._pixel_idx[int(gid)]
            Rmean = R_all[px].mean(axis=0)
            # Orthogonalise via SVD
            U, _, Vt = np.linalg.svd(Rmean)
            Rorth = U @ Vt
            if np.linalg.det(Rorth) < 0:
                U[:, -1] *= -1
                Rorth = U @ Vt
            # Extract Bunge angles
            PHI_m = np.arccos(np.clip(Rorth[2, 2], -1, 1))
            if abs(np.sin(PHI_m)) < 1e-8:
                phi1_m = np.arctan2(-Rorth[0, 1], Rorth[0, 0]) if Rorth[2, 2] > 0 else \
                          np.arctan2(Rorth[0, 1], -Rorth[0, 0])
                phi2_m = 0.0
            else:
                phi1_m = np.arctan2(Rorth[2, 0] / np.sin(PHI_m),
                                    -Rorth[2, 1] / np.sin(PHI_m))
                phi2_m = np.arctan2(Rorth[0, 2] / np.sin(PHI_m),
                                    Rorth[1, 2] / np.sin(PHI_m))
            out[idx] = np.degrees([phi1_m % (2 * np.pi),
                                    PHI_m,
                                    phi2_m % (2 * np.pi)])
        return out

    @cached_property
    def gos(self) -> np.ndarray:
        """
        Grain Orientation Spread (degrees):
        mean misorientation of each pixel from its grain's mean orientation.
        """
        phi1_all = np.radians(self.data['Euler1'].values)
        PHI_all  = np.radians(self.data['Euler2'].values)
        phi2_all = np.radians(self.data['Euler3'].values)
        R_all = euler_to_rotation_matrix(phi1_all, PHI_all, phi2_all)

        phi1_m = np.radians(self.mean_euler[:, 0])
        PHI_m  = np.radians(self.mean_euler[:, 1])
        phi2_m = np.radians(self.mean_euler[:, 2])
        R_mean = euler_to_rotation_matrix(phi1_m, PHI_m, phi2_m)

        out = np.zeros(self.n_grains, dtype=np.float64)
        for idx, gid in enumerate(self._valid_ids):
            px = self._pixel_idx[int(gid)]
            Rg = R_mean[idx]                         # (3,3) mean
            Rpx = R_all[px]                          # (K,3,3) pixels
            M = np.einsum('ij,kjl->kil', Rg.T, Rpx)  # Rg^T @ Rpx[k]
            trace = np.einsum('kii->k', M)
            angles = np.degrees(np.arccos(np.clip((trace - 1) / 2, -1, 1)))
            out[idx] = angles.mean()
        return out

    # ------------------------------------------------------------------
    # Neighbours
    # ------------------------------------------------------------------

    def neighbors(self) -> dict[int, list[int]]:
        """
        Return dict mapping grain_id → list of neighbouring grain IDs.
        Two grains are neighbours if they share at least one pixel edge.
        """
        if self._all_pairs is None:
            return {int(g): [] for g in self._valid_ids}

        g = self.grain_ids
        pairs = self._all_pairs
        g1 = g[pairs[:, 0]]
        g2 = g[pairs[:, 1]]
        diff = (g1 != g2) & (g1 >= 0) & (g2 >= 0)

        neigh: dict[int, set] = {int(gid): set() for gid in self._valid_ids}
        for a, b in zip(g1[diff], g2[diff]):
            neigh[int(a)].add(int(b))
            neigh[int(b)].add(int(a))
        return {k: sorted(v) for k, v in neigh.items()}

    # ------------------------------------------------------------------
    # Merge grains
    # ------------------------------------------------------------------

    def merge(self, threshold_deg: float = 5.0) -> 'Grain2D':
        """
        Merge neighbouring grains whose mean-orientation misorientation is
        below threshold_deg.  Returns a new Grain2D.
        """
        from scipy.sparse.csgraph import connected_components
        from scipy.sparse import csr_matrix

        phi1 = np.radians(self.mean_euler[:, 0])
        PHI  = np.radians(self.mean_euler[:, 1])
        phi2 = np.radians(self.mean_euler[:, 2])
        R_mean = euler_to_rotation_matrix(phi1, PHI, phi2)

        neigh = self.neighbors()
        # Build merge-graph
        id_to_idx = {int(g): i for i, g in enumerate(self._valid_ids)}
        rows, cols = [], []
        for gid, nlist in neigh.items():
            i = id_to_idx[gid]
            for nid in nlist:
                j = id_to_idx[nid]
                Ri, Rj = R_mean[i], R_mean[j]
                M = Ri.T @ Rj
                trace = np.clip(np.trace(M), -1, 3)
                angle_deg = np.degrees(np.arccos((trace - 1) / 2))
                if angle_deg < threshold_deg:
                    rows.append(i); cols.append(j)

        n = self.n_grains
        if rows:
            A = csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(n, n))
            A = A + A.T
            _, comp_labels = connected_components(A, directed=False)
        else:
            comp_labels = np.arange(n)

        # Build new grain_ids mapping
        new_ids = self.grain_ids.copy()
        for px_idx in range(len(self.data)):
            old = self.grain_ids[px_idx]
            if old < 0:
                continue
            grain_idx = id_to_idx[old]
            new_ids[px_idx] = comp_labels[grain_idx]

        # Renumber
        valid = np.unique(new_ids[new_ids >= 0])
        remap = np.full(int(new_ids.max()) + 2, -1, dtype=np.int32)
        remap[valid] = np.arange(len(valid))
        new_ids = np.where(new_ids >= 0, remap[new_ids], -1)

        return Grain2D(self.data, new_ids, self._all_pairs)

    # ------------------------------------------------------------------
    # Selection
    # ------------------------------------------------------------------

    def select_by_phase(self, phase_id: int) -> 'Grain2D':
        keep = self._valid_ids[self.phase_id == phase_id]
        new_ids = np.where(np.isin(self.grain_ids, keep), self.grain_ids, -1)
        return Grain2D(self.data, new_ids, self._all_pairs)

    def select_by_size(self, min_px: int = 1, max_px: int | None = None) -> 'Grain2D':
        mask = self.n_pixels >= min_px
        if max_px is not None:
            mask &= self.n_pixels <= max_px
        keep = self._valid_ids[mask]
        new_ids = np.where(np.isin(self.grain_ids, keep), self.grain_ids, -1)
        return Grain2D(self.data, new_ids, self._all_pairs)

    def select_by_orientation(self, ref_euler: np.ndarray,
                               tolerance_deg: float = 15.0) -> 'Grain2D':
        """Keep grains whose mean orientation is within tolerance of ref_euler (degrees)."""
        phi1_m = np.radians(self.mean_euler[:, 0])
        PHI_m  = np.radians(self.mean_euler[:, 1])
        phi2_m = np.radians(self.mean_euler[:, 2])
        Rg = euler_to_rotation_matrix(phi1_m, PHI_m, phi2_m)

        ref_r = np.radians(ref_euler)
        Rref = euler_to_rotation_matrix(ref_r[0:1], ref_r[1:2], ref_r[2:3])[0]

        angles = []
        for i in range(self.n_grains):
            M = Rref.T @ Rg[i]
            trace = np.clip(np.trace(M), -1, 3)
            angles.append(np.degrees(np.arccos((trace - 1) / 2)))
        angles = np.array(angles)

        keep = self._valid_ids[angles <= tolerance_deg]
        new_ids = np.where(np.isin(self.grain_ids, keep), self.grain_ids, -1)
        return Grain2D(self.data, new_ids, self._all_pairs)

    # ------------------------------------------------------------------
    # Grain boundary object
    # ------------------------------------------------------------------

    @cached_property
    def boundary(self) -> 'GrainBoundary':
        from .grain_boundary import GrainBoundary
        return GrainBoundary(self)

    # ------------------------------------------------------------------
    # Triple points
    # ------------------------------------------------------------------

    @cached_property
    def triple_points(self) -> 'TriplePoints':
        from .triple_points import TriplePoints
        return TriplePoints(self)

    # ------------------------------------------------------------------
    # Summary DataFrame
    # ------------------------------------------------------------------

    def summary_df(self) -> pd.DataFrame:
        ell = self.fit_ellipse()
        ch  = self.convex_hull_props()
        fer = self.feret_diameter()
        return pd.DataFrame({
            "grain_id":           self._valid_ids,
            "phase":              self.phase_id,
            "n_pixels":           self.n_pixels,
            "area":               self.area,
            "perimeter":          self.perimeter,
            "equiv_radius":       self.equivalent_radius,
            "equiv_perimeter":    self.equivalent_perimeter,
            "diameter":           self.diameter,
            "aspect_ratio":       self.aspect_ratio,
            "gos_deg":            self.gos,
            "ellipse_a":          ell["a"],
            "ellipse_b":          ell["b"],
            "ellipse_angle_deg":  ell["angle"],
            "centroid_x":         self.centroid[:, 0],
            "centroid_y":         self.centroid[:, 1],
            "hull_area":          ch["hull_area"],
            "hull_perimeter":     ch["hull_perimeter"],
            "solidity":           ch["solidity"],
            "max_feret":          fer["max_feret"],
            "min_feret":          fer["min_feret"],
            "euler1_mean_deg":    self.mean_euler[:, 0],
            "euler2_mean_deg":    self.mean_euler[:, 1],
            "euler3_mean_deg":    self.mean_euler[:, 2],
        })

    # ------------------------------------------------------------------
    # Quick plots (return axes so callers can customise)
    # ------------------------------------------------------------------

    def plot_phase_map(self, ax=None, cmap: str = "tab10"):
        import matplotlib.pyplot as plt
        from matplotlib.colors import BoundaryNorm
        if ax is None:
            _, ax = plt.subplots(figsize=(7, 5))
        phase_map = np.full(len(self.data), -1, dtype=np.int32)
        phase_map[:] = self.data['Phase'].values
        sc = ax.scatter(self.data['X'], self.data['Y'],
                        c=self.grain_ids, cmap=cmap, s=1, rasterized=True)
        ax.set_aspect('equal')
        ax.set_title("Grain map (colour = grain ID)")
        ax.set_xlabel("X"); ax.set_ylabel("Y")
        return ax

    def plot_grain_size_map(self, ax=None, cmap: str = "viridis"):
        import matplotlib.pyplot as plt
        if ax is None:
            _, ax = plt.subplots(figsize=(7, 5))
        size_map = np.full(len(self.data), np.nan)
        area_dict = {int(gid): self.area[i]
                     for i, gid in enumerate(self._valid_ids)}
        for idx, row_grain in enumerate(self.grain_ids):
            if row_grain >= 0:
                size_map[idx] = area_dict.get(int(row_grain), np.nan)
        sc = ax.scatter(self.data['X'], self.data['Y'],
                        c=size_map, cmap=cmap, s=1, rasterized=True)
        plt.colorbar(sc, ax=ax, label="Grain area (map units²)")
        ax.set_aspect('equal')
        ax.set_title("Grain size map")
        ax.set_xlabel("X"); ax.set_ylabel("Y")
        return ax


# ---------------------------------------------------------------------------
# Convenience constructor
# ---------------------------------------------------------------------------

def from_ebsd(data: pd.DataFrame,
              threshold_deg: float = 10.0,
              min_pixels: int = 2) -> Grain2D:
    """Reconstruct grains from an EBSD DataFrame and return a Grain2D object."""
    grain_ids, pairs = calc_grains(data, threshold_deg, min_pixels)
    return Grain2D(data, grain_ids, pairs)
