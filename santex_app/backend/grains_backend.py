"""
GrainsBackend — thin wrapper around santex.grains for the PyQt5 app.

Holds the current Grain2D state and exposes convenience methods for the UI.
"""

from __future__ import annotations
import numpy as np
import pandas as pd


class GrainsBackend:
    def __init__(self):
        self._grains = None      # santex.grains.Grain2D | None
        self._data: pd.DataFrame | None = None

    # ------------------------------------------------------------------
    # State
    # ------------------------------------------------------------------

    @property
    def is_ready(self) -> bool:
        return self._grains is not None

    @property
    def n_grains(self) -> int:
        return self._grains.n_grains if self._grains else 0

    # ------------------------------------------------------------------
    # Reconstruction
    # ------------------------------------------------------------------

    def reconstruct(self, data: pd.DataFrame,
                    threshold_deg: float = 10.0,
                    min_pixels: int = 2) -> None:
        """Run grain reconstruction on an EBSD DataFrame."""
        from santex.grains import from_ebsd
        self._data = data
        self._grains = from_ebsd(data, threshold_deg, min_pixels)

    def get_grain_ids(self) -> np.ndarray | None:
        return self._grains.grain_ids if self._grains else None

    # ------------------------------------------------------------------
    # Shape parameters
    # ------------------------------------------------------------------

    def summary_df(self) -> pd.DataFrame | None:
        return self._grains.summary_df() if self._grains else None

    def fit_ellipse(self) -> dict | None:
        return self._grains.fit_ellipse() if self._grains else None

    def convex_hull_props(self) -> dict | None:
        return self._grains.convex_hull_props() if self._grains else None

    def feret_diameter(self) -> dict | None:
        return self._grains.feret_diameter() if self._grains else None

    # ------------------------------------------------------------------
    # Orientation analysis
    # ------------------------------------------------------------------

    def calc_kam(self, max_angle_deg: float = 5.0) -> np.ndarray | None:
        if not self.is_ready:
            return None
        from santex.grains import calc_kam
        return calc_kam(self._grains, max_angle_deg=max_angle_deg)

    def calc_grod(self) -> np.ndarray | None:
        if not self.is_ready:
            return None
        from santex.grains import calc_grod
        return calc_grod(self._grains)

    def gos(self) -> np.ndarray | None:
        return self._grains.gos if self._grains else None

    def mean_euler(self) -> np.ndarray | None:
        return self._grains.mean_euler if self._grains else None

    def texture_index(self) -> float:
        if not self.is_ready:
            return 1.0
        from santex.grains import calc_texture_index
        return calc_texture_index(self._grains.mean_euler)

    # ------------------------------------------------------------------
    # Grain boundaries
    # ------------------------------------------------------------------

    def boundary(self):
        return self._grains.boundary if self._grains else None

    def boundary_summary(self) -> pd.DataFrame | None:
        if not self.is_ready:
            return None
        gb = self._grains.boundary
        return pd.DataFrame({
            "misorientation_angle": gb.misorientation_angle,
            "segment_length":       gb.segment_length,
            "grain_a":              gb.grain_id_pairs[:, 0],
            "grain_b":              gb.grain_id_pairs[:, 1],
        })

    def subgrain_boundary_fraction(self, lagb_threshold: float = 15.0) -> np.ndarray | None:
        if not self.is_ready:
            return None
        from santex.grains import subgrain_boundary_fraction
        return subgrain_boundary_fraction(self._grains, lagb_threshold)

    def detect_csl(self, threshold_deg: float = 5.0) -> np.ndarray | None:
        if not self.is_ready:
            return None
        return self._grains.boundary.sigma_value(threshold_deg)

    def detect_twins(self, threshold_deg: float = 5.0) -> np.ndarray | None:
        if not self.is_ready:
            return None
        return self._grains.boundary.is_twinning(3, threshold_deg)

    def misorientation_distribution(self, bins: int = 36) -> tuple | None:
        if not self.is_ready:
            return None
        from santex.grains import calc_misorientation_distribution
        return calc_misorientation_distribution(self._grains.boundary, bins)

    # ------------------------------------------------------------------
    # Triple points
    # ------------------------------------------------------------------

    def triple_points_df(self) -> pd.DataFrame | None:
        if not self.is_ready:
            return None
        return self._grains.triple_points.to_dataframe()

    # ------------------------------------------------------------------
    # Merge / filter
    # ------------------------------------------------------------------

    def merge(self, threshold_deg: float = 5.0) -> None:
        if self.is_ready:
            self._grains = self._grains.merge(threshold_deg)

    def select_by_phase(self, phase_id: int) -> None:
        if self.is_ready:
            self._grains = self._grains.select_by_phase(phase_id)

    def select_by_size(self, min_px: int, max_px: int | None = None) -> None:
        if self.is_ready:
            self._grains = self._grains.select_by_size(min_px, max_px)

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def export_csv(self, filepath: str) -> None:
        if self.is_ready:
            from santex.grains import to_csv
            to_csv(self._grains, filepath)

    def export_neper_ori(self, filepath: str) -> None:
        if self.is_ready:
            from santex.grains import to_neper_ori
            to_neper_ori(self._grains, filepath)

    def export_neper_tess(self, filepath: str) -> None:
        if self.is_ready:
            from santex.grains import to_neper_tess
            to_neper_tess(self._grains, filepath)

    def export_ctf(self, filepath: str) -> None:
        if self.is_ready:
            from santex.grains import to_ctf
            to_ctf(self._grains, filepath)
