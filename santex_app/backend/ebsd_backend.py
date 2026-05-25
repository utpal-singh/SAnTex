import numpy as np
import pandas as pd
from santex.ebsd import EBSD


# MTEX-style colours assigned in order to non-zero phases
_MTEX_COLORS = [
    "LightSkyBlue", "DarkSeaGreen", "Goldenrod", "LightCoral",
    "LightGreen", "Plum", "PeachPuff", "Thistle",
]

# Laue group number → symmetry label
_LAUE_LABEL = {
    1:  "-1",      # triclinic
    2:  "2/m",     # monoclinic
    3:  "mmm",     # orthorhombic
    5:  "4/mmm",   # tetragonal
    7:  "-3m",     # trigonal
    9:  "6/mmm",   # hexagonal
    10: "m-3",     # cubic Th
    11: "m-3m",    # cubic Oh
}

# Crystal reference frames for low-symmetry groups
_CRYSTAL_REF = {
    2:  "X||a*, Z||c*",
    7:  "X||a, Y||b, Z||c*",
    9:  "X||a, Y||b, Z||c",
}


class EBSDBackend:
    def __init__(self):
        self.ebsd: EBSD | None = None
        self.data: pd.DataFrame | None = None
        self._data_original: pd.DataFrame | None = None   # kept for ROI reset
        self.filename: str | None = None
        self._header_data: pd.DataFrame | None = None
        self._phases_df: pd.DataFrame | None = None

    # ------------------------------------------------------------------
    # File loading
    # ------------------------------------------------------------------

    def load(self, filename: str) -> None:
        self.filename = filename
        self.ebsd = EBSD(filename)
        self.data, self._header_data = self.ebsd.get_ebsd_data_header()
        self._data_original = self.data.copy()
        # Parse phases table from CTF
        try:
            from santex.ebsd._ctf_parser import Ctf
            ctf = Ctf(filename)
            _, _, _, _, _, _, self._phases_df = ctf.phases_info()
        except Exception:
            self._phases_df = None

    @property
    def is_loaded(self) -> bool:
        return self.ebsd is not None and self.data is not None

    # ------------------------------------------------------------------
    # Phase info
    # ------------------------------------------------------------------

    def phase_rows(self) -> list[tuple[int, str, float]]:
        """Return list of (index, name, percentage) tuples for *indexed* phases.

        Phase 0 (unindexed / NaN pixels) is intentionally excluded: those
        pixels carry no valid orientation and must never contribute to any
        elastic-tensor average.
        """
        if not self.is_loaded:
            return []
        return [(idx, name, pct)
                for idx, name, pct in self.ebsd.phases()
                if idx != 0]   # skip unindexed pixels

    def phase_names(self) -> list[str]:
        if not self.is_loaded:
            return []
        return list(self.ebsd.phases_names()["phase"])

    def phases_df(self) -> pd.DataFrame | None:
        if not self.is_loaded:
            return None
        return self.ebsd.phases_df()

    # ------------------------------------------------------------------
    # Euler angles
    # ------------------------------------------------------------------

    def euler_angles(self, phase_index: int) -> pd.DataFrame | None:
        if not self.is_loaded:
            return None
        return self.ebsd.get_euler_angles(phase_index, data=self.data)

    # ------------------------------------------------------------------
    # Map data
    # ------------------------------------------------------------------

    def map_data(self) -> pd.DataFrame | None:
        """Return X, Y, Phase columns for map rendering."""
        if self.data is None:
            return None
        return self.data[["X", "Y", "Phase"]].copy()

    # ------------------------------------------------------------------
    # Filtering
    # ------------------------------------------------------------------

    def filter_mad(self, threshold: float = 0.7) -> None:
        if self.data is not None:
            self.data = self.ebsd.filter_MAD(self.data, threshold)

    def downsample(self, factor: int = 10) -> pd.DataFrame | None:
        if self.data is None:
            return None
        return self.ebsd.downsample_EBSD(factor, self.data)

    # ------------------------------------------------------------------
    # Tensor averaging from EBSD orientations
    # ------------------------------------------------------------------

    def voigt_average(self, cij_list: list, density_list: list,
                      phase_indices: list[int]) -> tuple[np.ndarray | None, float | None]:
        """Simple Voigt average (arithmetic mean of rotated tensors)."""
        if not self.is_loaded:
            return None, None
        euler_list = [self.euler_angles(i) for i in phase_indices]
        return self.ebsd.get_anisotropy_for_ebsd(
            cij=cij_list, euler_angles=euler_list, density=density_list, method=0
        )

    def vrh_average(self, cij_list: list, density_list: list,
                    phase_indices: list[int],
                    method: str = "voigt") -> tuple[np.ndarray | None, float | None]:
        """Voigt / Reuss / Hill average."""
        if not self.is_loaded:
            return None, None
        euler_list = [self.euler_angles(i) for i in phase_indices]
        return self.ebsd.get_voigt_reuss_hill(
            cij=cij_list, euler_angles=euler_list,
            density=density_list, method=method
        )

    # ------------------------------------------------------------------
    # Pole / inverse-pole figures (delegates to crystal.py)
    # ------------------------------------------------------------------

    def pole_figure_data(self, phase_index: int, uvw: list[int] = [1, 0, 0],
                          crystal_symmetry: str = "D2h"):
        """Return orix vectors for a pole figure."""
        if not self.is_loaded:
            return None
        from santex.ebsd.crystal import pf
        return pf(self.data, phase=phase_index, uvw=uvw,
                  crystal_symmetry=crystal_symmetry, random_val=False)

    # ------------------------------------------------------------------
    # Scan geometry helpers
    # ------------------------------------------------------------------

    def get_scan_info(self) -> dict:
        """Return dict with XStep, YStep, XCells, YCells from CTF header."""
        info = {"XStep": 1.0, "YStep": 1.0, "XCells": 0, "YCells": 0}
        if self._header_data is not None:
            hd = self._header_data.set_index("info")["value"]
            for k in info:
                if k in hd.index:
                    try:
                        info[k] = float(hd[k])
                    except (ValueError, TypeError):
                        pass
        return info

    # ------------------------------------------------------------------
    # MTEX-style summary strings
    # ------------------------------------------------------------------

    def mtex_ebsd_summary(self) -> str:
        """Return a multi-line string mimicking MTEX's `ebsd` display."""
        if not self.is_loaded or self.data is None:
            return "No EBSD data loaded."

        lines = []
        lines.append("ebsd = EBSD  (y↑→x)\n")

        # ---- phase table ----
        total = len(self.data)
        col_width = [7, 18, 16, 12, 10, 26]
        hdr = (f"{'Phase':>{col_width[0]}}  {'Orientations':<{col_width[1]}}"
               f"{'Mineral':<{col_width[2]}}{'Color':<{col_width[3]}}"
               f"{'Symmetry':<{col_width[4]}}{'Crystal reference frame':<{col_width[5]}}")
        lines.append(hdr)
        lines.append("-" * sum(col_width))

        phase_ids = sorted(self.data["Phase"].unique())
        color_idx = 0
        for pid in phase_ids:
            count = int((self.data["Phase"] == pid).sum())
            pct = 100.0 * count / total
            if pid == 0:
                mineral = "notIndexed"
                color = ""
                sym = ""
                ref = ""
            else:
                if self._phases_df is not None and pid in self._phases_df.index:
                    row = self._phases_df.loc[pid]
                    mineral = str(row.get("phase", f"Phase{pid}")).strip()
                    try:
                        laue_num = int(float(str(row.get("symmetry", 11)).split()[0]))
                    except (ValueError, TypeError):
                        laue_num = 11
                    sym = _LAUE_LABEL.get(laue_num, str(laue_num))
                    ref = _CRYSTAL_REF.get(laue_num, "")
                else:
                    mineral = f"Phase{pid}"
                    sym = ""
                    ref = ""
                color = _MTEX_COLORS[color_idx % len(_MTEX_COLORS)]
                color_idx += 1

            lines.append(
                f"{pid:>{col_width[0]}}  {f'{count} ({pct:.1f}%)':<{col_width[1]}}"
                f"{mineral:<{col_width[2]}}{color:<{col_width[3]}}"
                f"{sym:<{col_width[4]}}{ref:<{col_width[5]}}"
            )

        lines.append("")

        # ---- properties ----
        prop_cols = [c for c in self.data.columns
                     if c not in ("Phase", "X", "Y", "Euler1", "Euler2", "Euler3")]
        if prop_cols:
            lines.append(f" Properties: {', '.join(c.lower() for c in prop_cols)}")

        # ---- scan info ----
        si = self.get_scan_info()
        lines.append(" Scan unit : um")

        x_vals = self.data["X"]
        y_vals = self.data["Y"]
        lines.append(f" X x Y x Z : [{x_vals.min():.0f}, {x_vals.max():.0f}]"
                     f" x [{y_vals.min():.0f}, {y_vals.max():.0f}] x [0, 0]")
        lines.append(" Normal vector: (0,0,1)")

        return "\n".join(lines)

    def mtex_grains_summary(self, grains_backend) -> str:
        """Return a multi-line string mimicking MTEX's `grain2d` display."""
        if not grains_backend.is_ready:
            return "No grains reconstructed yet."

        g = grains_backend._grains
        grain_ids_arr = g.grain_ids          # per-pixel grain id array
        lines = []
        lines.append("grains = grain2d  (y↑→x)\n")

        # ---- phase table ----
        col_width = [7, 8, 8, 16, 12, 26]
        hdr = (f"{'Phase':>{col_width[0]}}  {'Grains':<{col_width[1]}}"
               f"{'Pixels':<{col_width[2]}}{'Mineral':<{col_width[3]}}"
               f"{'Symmetry':<{col_width[4]}}{'Crystal reference frame':<{col_width[5]}}")
        lines.append(hdr)
        lines.append("-" * sum(col_width))

        data = grains_backend._data if grains_backend._data is not None else self.data
        phase_ids = sorted(data["Phase"].unique()) if data is not None else []
        color_idx = 0
        for pid in phase_ids:
            phase_mask = (data["Phase"] == pid) if data is not None else np.array([], dtype=bool)
            px_count = int(phase_mask.sum())
            # count unique positive grain ids in this phase
            phase_grain_ids = grain_ids_arr[phase_mask.values] if data is not None else np.array([])
            unique_grains = len(set(phase_grain_ids[phase_grain_ids >= 0]))

            if pid == 0:
                mineral = "notIndexed"
                sym = ""
                ref = ""
            else:
                if self._phases_df is not None and pid in self._phases_df.index:
                    row = self._phases_df.loc[pid]
                    mineral = str(row.get("phase", f"Phase{pid}")).strip()
                    try:
                        laue_num = int(float(str(row.get("symmetry", 11)).split()[0]))
                    except (ValueError, TypeError):
                        laue_num = 11
                    sym = _LAUE_LABEL.get(laue_num, str(laue_num))
                    ref = _CRYSTAL_REF.get(laue_num, "")
                else:
                    mineral = f"Phase{pid}"
                    sym = ""
                    ref = ""
                color_idx += 1

            lines.append(
                f"{pid:>{col_width[0]}}  {unique_grains:<{col_width[1]}}"
                f"{px_count:<{col_width[2]}}{mineral:<{col_width[3]}}"
                f"{sym:<{col_width[4]}}{ref:<{col_width[5]}}"
            )

        lines.append("")

        # ---- boundary stats ----
        try:
            si = self.get_scan_info()
            step = si.get("XStep", 1.0)
            gb = g.boundary
            n_seg = len(gb.misorientation_angle)
            total_um = n_seg * step
            lines.append(f" boundary segments: {n_seg}  ({total_um / 1e3:.2g}e+03 µm)"
                         if total_um >= 1000 else
                         f" boundary segments: {n_seg}  ({total_um:.0f} µm)")
        except Exception:
            lines.append(" boundary segments: N/A")

        try:
            inner = int((g.boundary.misorientation_angle == 0).sum())
            lines.append(f" inner boundary segments: {inner}")
        except Exception:
            pass

        try:
            tp = g.triple_points
            lines.append(f" triple points: {len(tp.x)}")
        except Exception:
            pass

        lines.append("")
        lines.append(" Properties: meanRotation, GOS")
        return "\n".join(lines)

    # ------------------------------------------------------------------
    # ROI clipping
    # ------------------------------------------------------------------

    def clip_roi(self, x_min: float, x_max: float,
                 y_min: float, y_max: float) -> None:
        """Clip data to a rectangular region of interest."""
        if self._data_original is None:
            return
        d = self._data_original
        mask = (
            (d["X"] >= x_min) & (d["X"] <= x_max) &
            (d["Y"] >= y_min) & (d["Y"] <= y_max)
        )
        self.data = d[mask].reset_index(drop=True)

    def reset_roi(self) -> None:
        """Restore the full dataset (undo ROI clip)."""
        if self._data_original is not None:
            self.data = self._data_original.copy()

    # ------------------------------------------------------------------

    @property
    def n_points(self) -> int:
        return len(self.data) if self.data is not None else 0
