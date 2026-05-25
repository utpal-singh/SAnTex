"""
ODFBackend — thin wrapper around santex.odf for the PyQt5 app.

Holds the current ODF state and exposes methods for the UI.
Crystal symmetry names are resolved through santex.grains.symmetry,
which supports all 32 point groups, all 11 Laue groups, Schoenflies
symbols, HM symbols, and MTEX-compatible aliases.
"""

from __future__ import annotations
import numpy as np
import pandas as pd


class ODFBackend:
    def __init__(self):
        self._odf = None           # santex.odf.ODF | None
        self._euler_deg = None     # (N, 3) cached Euler angles
        self._phase_id = None      # which phase was used
        self.crystal_symmetry = "orthorhombic"
        self.halfwidth_deg = 10.0
        # cached results
        self._sections: dict | None = None
        self._pf_vectors: dict | None = None  # hkl_str → (vectors, weights)
        self._ipf_colors: np.ndarray | None = None

    # ------------------------------------------------------------------
    # State
    # ------------------------------------------------------------------

    @property
    def is_ready(self) -> bool:
        return self._odf is not None

    # ------------------------------------------------------------------
    # ODF Computation
    # ------------------------------------------------------------------

    @staticmethod
    def symmetry_display_list() -> list[tuple[str, str]]:
        """
        Return a list of (display_label, canonical_key) tuples for the UI dropdown.

        Ordered by Laue group id (triclinic → cubic).
        Returns tuples like ("Cubic (m-3m)  —  24 ops", "m-3m").
        """
        from santex.grains.symmetry import list_laue_groups
        result = []
        for lg in list_laue_groups():
            label = (f"{lg['display']}"
                     f"  [{lg['schoenflies']}]"
                     f"  —  {lg['proper_order']} ops")
            result.append((label, lg["laue_group"]))
        return result

    def compute_odf(self,
                    data: pd.DataFrame,
                    phase: int = 1,
                    crystal_symmetry: str = "mmm",
                    halfwidth_deg: float = 10.0,
                    max_orientations: int = 5000) -> None:
        """
        Compute ODF from EBSD data for a given phase.

        Parameters
        ----------
        data              : pd.DataFrame with Phase, Euler1, Euler2, Euler3
        phase             : EBSD phase index to use (0 = unindexed, skip)
        crystal_symmetry  : any recognised symmetry name/alias
                            (e.g. 'mmm', 'orthorhombic', 'm-3m', 'cubic',
                             'Oh', '-3m', 'trigonal', '2/m', …)
        halfwidth_deg     : kernel halfwidth in degrees
        max_orientations  : subsample cap
        """
        from santex.odf import odf_from_ebsd
        from santex.grains.symmetry import get_laue_group

        # Normalise to canonical Laue-group HM symbol (e.g. "mmm", "m-3m")
        laue_key = get_laue_group(crystal_symmetry)
        self.crystal_symmetry = laue_key
        self.halfwidth_deg    = halfwidth_deg
        self._phase_id        = phase

        self._odf = odf_from_ebsd(
            data, phase=phase,
            crystal_symmetry=laue_key,
            halfwidth_deg=halfwidth_deg,
            max_orientations=max_orientations,
        )
        # Cache Euler angles for other computations
        subset = data[data["Phase"] == phase] if phase is not None else data[data["Phase"] > 0]
        self._euler_deg = subset[["Euler1", "Euler2", "Euler3"]].values

        # Invalidate cached results
        self._sections   = None
        self._pf_vectors = None
        self._ipf_colors = None

    # ------------------------------------------------------------------
    # ODF Sections
    # ------------------------------------------------------------------

    def get_phi2_sections(self,
                          phi2_values_deg: list[float] | None = None,
                          resolution_deg: float = 5.0
                          ) -> dict | None:
        """Return phi2-section dict: phi2 → (phi1_2d, PHI_2d, odf_2d)."""
        if not self.is_ready:
            return None
        if phi2_values_deg is None:
            phi2_values_deg = [0, 15, 30, 45, 60, 75, 90]
        return self._odf.phi2_sections(phi2_values_deg, resolution_deg)

    def get_sigma_sections(self,
                           sigma_values_deg: list[float] | None = None,
                           resolution_deg: float = 5.0) -> dict | None:
        if not self.is_ready:
            return None
        if sigma_values_deg is None:
            sigma_values_deg = [0, 15, 30, 45, 60, 75, 90]
        return self._odf.sigma_sections(sigma_values_deg, resolution_deg)

    def get_phi1_sections(self,
                          phi1_values_deg: list[float] | None = None,
                          resolution_deg: float = 5.0) -> dict | None:
        if not self.is_ready:
            return None
        if phi1_values_deg is None:
            phi1_values_deg = [0, 15, 30, 45, 60, 75, 90]
        return self._odf.phi1_sections(phi1_values_deg, resolution_deg)

    # ------------------------------------------------------------------
    # Pole Figures
    # ------------------------------------------------------------------

    def get_pole_figure(self,
                        hkl: tuple[int, int, int] = (0, 0, 1),
                        sample_symmetry: str = "none") -> tuple | None:
        """
        Compute pole figure vectors for a crystal direction.

        Returns (vectors (N,3), weights or None).
        """
        if self._euler_deg is None:
            return None
        from santex.odf import calc_pole_figure
        return calc_pole_figure(
            self._euler_deg, hkl,
            self.crystal_symmetry, sample_symmetry
        )

    def get_pole_figure_kde(self,
                            hkl: tuple = (0, 0, 1),
                            sample_symmetry: str = "none",
                            bandwidth_deg: float = 7.5
                            ) -> tuple | None:
        """
        Compute smooth pole figure (KDE on S²).

        Returns (xy_grid (M,2), density (M,), xy_scatter (N,2)).
        """
        if self._euler_deg is None:
            return None
        from santex.odf import calc_pole_figure, pole_figure_grid, sphere_kde, stereo_project
        vectors, weights = calc_pole_figure(
            self._euler_deg, hkl, self.crystal_symmetry, sample_symmetry
        )
        # Upper hemisphere only
        vectors_up = vectors.copy()
        neg = vectors_up[:, 2] < 0
        vectors_up[neg] = -vectors_up[neg]

        xy_grid, eval_vecs = pole_figure_grid(n_pts=2000)
        density = sphere_kde(vectors_up, eval_vecs,
                             bandwidth=np.radians(bandwidth_deg), weights=weights)
        xy_scatter = stereo_project(vectors_up)
        return xy_grid, density, xy_scatter

    # ------------------------------------------------------------------
    # Inverse Pole Figures
    # ------------------------------------------------------------------

    def get_ipf_colors(self,
                       sample_dir: tuple = (0, 0, 1)) -> np.ndarray | None:
        """
        Compute per-point IPF colours for the stored Euler angles.

        Returns (N, 3) RGB array in [0, 1].
        """
        if self._euler_deg is None:
            return None
        from santex.odf import ipf_map_colors
        return ipf_map_colors(self._euler_deg, sample_dir, self.crystal_symmetry)

    def get_ipf_triangle_colors(self,
                                n_pts: int = 5000) -> tuple[np.ndarray, np.ndarray] | None:
        """
        Generate IPF colour triangle for display.

        Returns (xy, rgb) where xy is in stereographic coords.
        """
        if not self.is_ready:
            return None
        from santex.odf import stereo_project, ipf_color
        from santex.odf.pole_figure import _into_fundamental_sector, _symmetry_ops

        csym = _symmetry_ops(self.crystal_symmetry)

        # Sample directions in the fundamental sector
        rng = np.random.default_rng(0)
        theta = rng.uniform(0, np.pi / 2, n_pts)
        phi   = rng.uniform(0, np.pi / 2, n_pts)
        x = np.sin(theta) * np.cos(phi)
        y = np.sin(theta) * np.sin(phi)
        z = np.cos(theta)
        dirs = np.stack([x, y, z], axis=1)
        # Bring into fundamental sector
        dirs = _into_fundamental_sector(dirs, csym, self.crystal_symmetry)

        rgb  = ipf_color(dirs, self.crystal_symmetry)
        xy   = stereo_project(dirs)
        return xy, rgb

    # ------------------------------------------------------------------
    # ODF Properties
    # ------------------------------------------------------------------

    def texture_index(self, resolution_deg: float = 10.0) -> float | None:
        if not self.is_ready:
            return None
        return self._odf.texture_index(resolution_deg)

    def entropy(self, resolution_deg: float = 10.0) -> float | None:
        if not self.is_ready:
            return None
        return self._odf.entropy(resolution_deg)

    def odf_max(self, resolution_deg: float = 5.0) -> float | None:
        if not self.is_ready:
            return None
        return self._odf.odf_max(resolution_deg)

    # ------------------------------------------------------------------
    # Texture Components
    # ------------------------------------------------------------------

    def component_analysis(self, max_angle_deg: float = 15.0) -> list[dict]:
        if not self.is_ready:
            return []
        return self._odf.component_analysis(max_angle_deg)

    # ------------------------------------------------------------------
    # Random sampling
    # ------------------------------------------------------------------

    def random_sampling(self, n: int = 1000) -> np.ndarray | None:
        if not self.is_ready:
            return None
        return self._odf.random_sampling(n)

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def export_popla(self, filepath: str, resolution_deg: float = 5.0) -> None:
        if self.is_ready:
            from santex.odf import to_popla
            to_popla(self._odf, filepath, resolution_deg)

    def export_mtex_txt(self, filepath: str, resolution_deg: float = 5.0) -> None:
        if self.is_ready:
            from santex.odf import to_mtex_txt
            to_mtex_txt(self._odf, filepath, resolution_deg)

    def export_vpsc(self, filepath: str, n: int = 1000) -> None:
        if self.is_ready:
            from santex.odf import to_vpsc
            to_vpsc(self._odf, filepath, n)
