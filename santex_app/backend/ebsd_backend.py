import numpy as np
import pandas as pd
from santex.ebsd import EBSD
from santex.ebsd.ebsd_operations import (
    ipf_map_colors, bc_bs_map_data,
    compute_kam, compute_mis2mean,
    denoise_mean_filter, fill_missing_data, regrid,
    select_by_condition, get_line_profile,
    orientation_scatter_data,
    simulate_ebsd,
    export_ctf, export_ang, export_csv, export_hdf5,
)
from santex.ebsd.ipf_coloring import spacegroup_to_laue


# MTEX-style colours assigned in order to non-zero phases
_MTEX_COLORS = [
    "LightSkyBlue", "DarkSeaGreen", "Goldenrod", "LightCoral",
    "LightGreen", "Plum", "PeachPuff", "Thistle",
]

# Laue-group key → Hermann–Mauguin symbol
_LAUE_KEY_LABEL = {
    "C1":  "1",
    "Ci":  "-1",
    "C2h": "2/m",
    "D2h": "mmm",
    "C4h": "4/m",
    "D4h": "4/mmm",
    "S6":  "-3",
    "D3d": "-3m",
    "C6h": "6/m",
    "D6h": "6/mmm",
    "Th":  "m-3",
    "Oh":  "m-3m",
}

# Laue-group key → crystal reference frame note (non-cubic / non-orthorhombic)
_LAUE_KEY_REF = {
    "Ci":  "X||a*, Z||c",
    "C2h": "X||a*, Z||c*",
    "S6":  "X||a, Z||c*",
    "D3d": "X||a, Z||c*",
    "C6h": "X||a, Z||c",
    "D6h": "X||a, Z||c",
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

    def load_from_project(self, project: "EBSDProject", dataset_id: str) -> None:
        """Load a dataset from an EBSDProject using its LRU cache.

        This avoids re-reading the CTF file from disk when the dataset is
        already cached.  The working copy (``self.data``) is separate from
        the cached raw DataFrame so that filters / ROI clips do not
        contaminate the cache entry for other backends.

        Parameters
        ----------
        project     : EBSDProject
        dataset_id  : str  UUID of the dataset to load
        """
        from backend.ebsd_project import EBSDProject   # local import avoids circularity

        meta = project.get_meta(dataset_id)
        if meta is None:
            raise ValueError(f"Dataset {dataset_id!r} not found in project.")

        # Build EBSD wrapper (fast — no file read in __init__)
        self.ebsd    = EBSD(meta.path)
        self.filename = meta.path

        # Get data from LRU cache (loads from disk only on cache miss)
        raw_df = project.get_data(dataset_id)
        if raw_df is None:
            raise RuntimeError(f"Could not load data for dataset {dataset_id!r}.")

        self.data             = raw_df.copy()   # working copy (mutable)
        self._data_original   = raw_df.copy()   # for ROI reset
        self._phases_df       = meta.phases_df
        self._header_data     = meta.header_data

    @property
    def is_loaded(self) -> bool:
        return self.ebsd is not None and self.data is not None

    # ------------------------------------------------------------------
    # Phase info
    # ------------------------------------------------------------------

    def phase_symmetry(self, phase_index: int) -> str:
        """Return the canonical Laue-group key for *phase_index* (auto-detected from CTF).

        Reads the ``symmetry`` field (ITA space-group number 1–230, or Channel-5
        Laue-group code 1–11) from the CTF phases table and converts it to a
        canonical Laue key such as ``'D2h'``, ``'Oh'``, ``'D6h'``, etc.

        Falls back to ``'D2h'`` (orthorhombic) when the phase is not found.
        """
        if self._phases_df is not None and phase_index in self._phases_df.index:
            sg = self._phases_df.loc[phase_index].get("symmetry", 11)
            return spacegroup_to_laue(sg)
        return "D2h"

    def phase_symmetries(self) -> dict[int, str]:
        """Return {phase_index: laue_key} for every indexed phase in the dataset."""
        if not self.is_loaded or self.data is None:
            return {}
        return {
            idx: self.phase_symmetry(idx)
            for idx in sorted(self.data["Phase"].unique())
            if idx != 0
        }

    def phase_symmetry_label(self, phase_index: int) -> str:
        """Human-readable symmetry string, e.g. ``'mmm (D2h)'``."""
        key = self.phase_symmetry(phase_index)
        hm  = _LAUE_KEY_LABEL.get(key, key)
        return f"{hm} ({key})"

    def phase_rows(self) -> list[tuple[int, str, float]]:
        """Return list of (index, name, percentage) tuples for *indexed* phases.

        Phase 0 (unindexed / NaN pixels) is intentionally excluded: those
        pixels carry no valid orientation and must never contribute to any
        elastic-tensor average.

        Passes ``self.data`` to ``EBSD.phases()`` so that the CTF file is
        never re-read from disk (important when data comes from the LRU cache).
        """
        if not self.is_loaded:
            return []
        return [(idx, name, pct)
                for idx, name, pct in self.ebsd.phases(df_phases=self.data)
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
                    laue_key = spacegroup_to_laue(row.get("symmetry", 11))
                    sym = _LAUE_KEY_LABEL.get(laue_key, laue_key)
                    ref = _LAUE_KEY_REF.get(laue_key, "")
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
                    laue_key = spacegroup_to_laue(row.get("symmetry", 11))
                    sym = _LAUE_KEY_LABEL.get(laue_key, laue_key)
                    ref = _LAUE_KEY_REF.get(laue_key, "")
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
    # Specimen reference frame correction  (MTEX-style)
    # ------------------------------------------------------------------

    def apply_reference_frame_correction(
        self,
        method: str,
        euler_rotation: tuple[float, float, float] | None = None,
    ) -> None:
        """Apply a specimen / Euler reference-frame correction in-place.

        Parameters
        ----------
        method : str
            ``'none'``            – no correction (identity)
            ``'euler2spatial'``   – add 180° to φ₂ (Euler3); converts the Euler
                                    reference frame to match the spatial (Y-down)
                                    frame used by Oxford/HKL CTF files.
                                    Equivalent to MTEX ``convertEuler2SpatialReferenceFrame``.
            ``'spatial2euler'``   – flip the Y spatial coordinate (y ← y_max − y);
                                    converts the spatial scan frame to match the
                                    Euler (Y-up) reference frame.
                                    Equivalent to MTEX ``convertSpatial2EulerReferenceFrame``.
            ``'custom_euler'``    – apply an arbitrary ZXZ Bunge rotation
                                    (given by *euler_rotation*) to every pixel's
                                    orientation, expressed as a post-multiplication
                                    R_new = R_pixel · R_custom.
            ``'custom_spatial'``  – rotate the spatial (X, Y) coordinates by the
                                    Z-angle taken from the first element of
                                    *euler_rotation* (the in-plane rotation).
        euler_rotation : (φ₁, Φ, φ₂) in **degrees**, required for ``custom_*``.
        """
        if self._data_original is None:
            return

        d = self._data_original

        if method == 'none':
            pass  # nothing to do

        elif method == 'euler2spatial':
            # φ₂ += 180° — rotates orientation by Rz(π) on the right (specimen side)
            d = d.copy()
            d['Euler3'] = (d['Euler3'] + 180.0) % 360.0

        elif method == 'spatial2euler':
            # Flip Y: y ← y_max − y
            d = d.copy()
            y_max = d['Y'].max()
            d['Y'] = y_max - d['Y']

        elif method in ('custom_euler', 'custom_spatial'):
            if euler_rotation is None:
                return
            phi1_c, Phi_c, phi2_c = [np.deg2rad(a) for a in euler_rotation]

            if method == 'custom_euler':
                # Build the custom rotation matrix (ZXZ Bunge)
                from santex.tensor.tensor import Tensor
                t = Tensor()
                R_custom = t.euler_to_rotation(phi1_c, Phi_c, phi2_c)

                # For each pixel: R_new = R_pixel · R_custom
                # Express R_new back as Bunge ZXZ angles.
                d = d.copy()
                new_e1, new_e, new_e2 = [], [], []
                for _, row in d.iterrows():
                    R_pix = t.euler_to_rotation(
                        np.deg2rad(row['Euler1']),
                        np.deg2rad(row['Euler2']),
                        np.deg2rad(row['Euler3']),
                    )
                    R_new = R_pix @ R_custom
                    e1, e, e2 = _rotation_matrix_to_bunge_zxz(R_new)
                    new_e1.append(np.rad2deg(e1))
                    new_e.append(np.rad2deg(e))
                    new_e2.append(np.rad2deg(e2))
                d['Euler1'] = new_e1
                d['Euler2'] = new_e
                d['Euler3'] = new_e2

            else:  # 'custom_spatial'
                # In-plane rotation only (around Z by φ₁)
                theta = phi1_c  # use first Euler angle as the in-plane rotation
                d = d.copy()
                cx, cy = d['X'].mean(), d['Y'].mean()
                dx, dy = d['X'] - cx, d['Y'] - cy
                cos_t, sin_t = np.cos(theta), np.sin(theta)
                d['X'] = cx + cos_t * dx - sin_t * dy
                d['Y'] = cy + sin_t * dx + cos_t * dy

        else:
            raise ValueError(f"Unknown reference-frame method: {method!r}")

        self._data_original = d
        self.data = d.copy()

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
    # IPF map colours
    # ------------------------------------------------------------------

    def ipf_map_colors(self, phase_index: int = 1,
                       direction: str = "ND",
                       color_key: str = "tsl") -> dict | None:
        """Per-pixel IPF RGB colours for *phase_index*.

        Crystal symmetry is auto-detected from the CTF phases table
        (``symmetry`` field → ITA space-group number → Laue group).

        Parameters
        ----------
        phase_index : int
            Phase number (> 0).
        direction   : ``'ND'`` (Z), ``'RD'`` (X), or ``'TD'`` (Y).
        color_key   : ``'tsl'`` (blue-green-red corner interpolation, MTEX default)
                      or ``'hsv'`` (white-at-pole, rainbow at boundary).
        """
        if self.data is None:
            return None
        sym = self.phase_symmetry(phase_index)
        return ipf_map_colors(self.data, phase_index, direction, sym, color_key)

    def bc_bs_map_data(self) -> dict | None:
        """Band-contrast / band-slope map data for all indexed pixels."""
        if self.data is None:
            return None
        return bc_bs_map_data(self.data)

    # ------------------------------------------------------------------
    # Orientation analysis (no grain reconstruction needed)
    # ------------------------------------------------------------------

    def compute_kam(self, kernel_size: int = 1,
                    max_angle_deg: float = 5.0,
                    same_phase_only: bool = True) -> np.ndarray | None:
        """Kernel Average Misorientation per pixel (degrees)."""
        if self.data is None:
            return None
        return compute_kam(self.data, kernel_size, max_angle_deg, same_phase_only)

    def compute_mis2mean(self, phase_index: int = 1) -> np.ndarray | None:
        """Per-pixel misorientation from the mean orientation of *phase_index*."""
        if self.data is None:
            return None
        return compute_mis2mean(self.data, phase_index)

    def orientation_scatter_data(self, phase_index: int = 1,
                                  axes: tuple = ("Euler1", "Euler2"),
                                  max_points: int = 5000) -> dict | None:
        """Euler-space scatter data for the given phase."""
        if self.data is None:
            return None
        return orientation_scatter_data(self.data, phase_index, axes, max_points)

    # ------------------------------------------------------------------
    # Pre-processing: denoising, fill, regrid
    # ------------------------------------------------------------------

    def denoise(self, kernel_size: int = 1,
                max_angle_deg: float = 5.0) -> None:
        """Apply mean-filter denoising to orientations (modifies ``data`` in-place)."""
        if self.data is None:
            return
        self.data = denoise_mean_filter(self.data, kernel_size, max_angle_deg)

    def fill_missing(self, method: str = "nearest") -> None:
        """Fill phase-0 pixels by nearest-neighbour / local-mean interpolation."""
        if self.data is None:
            return
        self.data = fill_missing_data(self.data, method)

    def regrid_data(self, target_step: float | None = None,
                    scale: float | None = None) -> None:
        """Resample the scan onto a new regular Cartesian grid."""
        if self.data is None:
            return
        self.data = regrid(self.data, target_step=target_step, scale=scale)

    # ------------------------------------------------------------------
    # Selection
    # ------------------------------------------------------------------

    def select_by_condition(self, **kwargs) -> pd.DataFrame | None:
        """Return a filtered copy of ``data``; does not modify the original."""
        if self.data is None:
            return None
        return select_by_condition(self.data, **kwargs)

    # ------------------------------------------------------------------
    # Line profile
    # ------------------------------------------------------------------

    def get_line_profile(self, x0: float, y0: float,
                         x1: float, y1: float,
                         n_points: int = 100,
                         scalars: list[str] | None = None) -> dict | None:
        """Extract a line profile through the EBSD map."""
        if self.data is None:
            return None
        return get_line_profile(self.data, x0, y0, x1, y1, n_points, scalars)

    # ------------------------------------------------------------------
    # Simulation
    # ------------------------------------------------------------------

    def simulate_from_euler(self, euler_angles: np.ndarray,
                             n_cols: int = 50, n_rows: int = 50,
                             step: float = 1.0, phase_id: int = 1,
                             noise_deg: float = 0.5) -> pd.DataFrame:
        """Generate a synthetic EBSD DataFrame from prototype orientations."""
        return simulate_ebsd(euler_angles, n_cols, n_rows, step, phase_id, noise_deg)

    def load_synthetic(self, df: pd.DataFrame) -> None:
        """Replace the current dataset with a synthetic one."""
        self.data = df.copy()
        self._data_original = df.copy()
        self.filename = "<synthetic>"

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def export(self, filepath: str, fmt: str = "csv") -> None:
        """Export the current EBSD data.

        Parameters
        ----------
        filepath : output file path
        fmt      : ``'csv'``, ``'ctf'``, ``'ang'``, or ``'hdf5'``
        """
        if self.data is None:
            return
        fmt = fmt.lower()
        if fmt == "ctf":
            export_ctf(self.data, filepath)
        elif fmt == "ang":
            export_ang(self.data, filepath)
        elif fmt == "hdf5":
            export_hdf5(self.data, filepath)
        else:
            export_csv(self.data, filepath)

    # ------------------------------------------------------------------

    @property
    def n_points(self) -> int:
        return len(self.data) if self.data is not None else 0


# ---------------------------------------------------------------------------
# Module-level helpers
# ---------------------------------------------------------------------------

def _rotation_matrix_to_bunge_zxz(R: np.ndarray) -> tuple[float, float, float]:
    """Convert a 3×3 rotation matrix to Bunge ZXZ Euler angles (radians).

    Returns (φ₁, Φ, φ₂) in [0, 2π) × [0, π] × [0, 2π).
    """
    # Φ from R[2,2] = cos(Φ)
    Phi = np.arccos(np.clip(R[2, 2], -1.0, 1.0))

    if np.abs(np.sin(Phi)) < 1e-10:
        # Gimbal lock: Φ ≈ 0 or π → φ₁ and φ₂ are not uniquely determined
        if Phi < 1.0:   # Φ ≈ 0
            phi1 = np.arctan2(R[0, 1], R[0, 0])
            phi2 = 0.0
        else:            # Φ ≈ π
            phi1 = np.arctan2(R[0, 1], -R[0, 0])
            phi2 = 0.0
    else:
        phi1 = np.arctan2(R[2, 0], -R[2, 1])
        phi2 = np.arctan2(R[0, 2],  R[1, 2])

    # Map to [0, 2π)
    phi1 %= (2 * np.pi)
    phi2 %= (2 * np.pi)
    return phi1, Phi, phi2
