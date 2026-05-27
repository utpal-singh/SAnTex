"""
EBSDTab — EBSD map visualisation with Plotly and per-phase colour pickers.

Plots
-----
  Phase map   : scatter coloured by phase (discrete, user-configurable colours)
  Pole figure : quick stereonet scatter from raw Euler data
  Stereonet   : anisotropy Vp/Vs after VRH averaging

ROI
---
  Spinbox-based X/Y range; use Plotly zoom + hover to read coordinates,
  then enter them in the spin boxes and click "Clip to ROI".
"""

from __future__ import annotations
import numpy as np

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QListWidget, QPushButton, QSplitter,
    QDoubleSpinBox, QSpinBox, QFormLayout, QComboBox,
    QProgressBar, QTableWidget, QTableWidgetItem, QHeaderView,
    QTextEdit, QTabWidget, QMessageBox, QSizePolicy, QColorDialog,
    QScrollArea, QCheckBox, QDialog, QDialogButtonBox, QRadioButton,
    QButtonGroup, QFrame,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QFont, QColor

from frontend.widgets.plotly_widget import (
    PlotlyWidget, COLORSCALES, DEFAULT_PHASE_COLORS, ColorscaleComboBox,
)
from frontend.tabs._stereonet import make_stereonet_figure, STEREONET_STYLES
from frontend.widgets.pyvista_widget import PyVistaWidget

import plotly.graph_objects as go


# ---------------------------------------------------------------------------
# Helper: convert RGBA uint8 numpy array → base-64 PNG data URI
# (used to embed the IPF colour-key legend as a Plotly layout image)
# ---------------------------------------------------------------------------

def _rgba_to_data_uri(rgba: "np.ndarray") -> str:
    """Encode a (H, W, 4) uint8 RGBA array as a ``data:image/png;base64,...`` URI."""
    import io, base64
    try:
        from PIL import Image as _PILImage
        buf = io.BytesIO()
        _PILImage.fromarray(rgba, mode="RGBA").save(buf, format="PNG")
        b64 = base64.b64encode(buf.getvalue()).decode()
        return f"data:image/png;base64,{b64}"
    except ImportError:
        # Pillow not available — return a tiny transparent placeholder
        _EMPTY_PNG = (
            "iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAYAAAAfFcSJAAAAC0lEQVQI12NgAAIA"
            "BQAABjE+ibYAAAAASUVORK5CYII="
        )
        return f"data:image/png;base64,{_EMPTY_PNG}"


# ---------------------------------------------------------------------------
# IPF figure builder  — map + labelled colour-key triangle (MTEX style)
# ---------------------------------------------------------------------------

def _build_ipf_figure(result: dict, phase_idx: int,
                      direction_str: str, color_key: str,
                      sym_label: str):
    """Build a Plotly figure with the IPF map on the left and the labelled
    colour-key triangle on the right, matching the MTEX IPF display style.

    Parameters
    ----------
    result       : dict returned by ``EBSDBackend.ipf_map_colors()``
                   keys: x, y, r, g, b, rgb_hex, sym_key
    phase_idx    : phase number (for title)
    direction_str: 'ND', 'RD', or 'TD'
    color_key    : 'tsl' or 'hsv'
    sym_label    : human-readable symmetry, e.g. 'mmm (D2h)'
    """
    from plotly.subplots import make_subplots
    from santex.ebsd.ipf_coloring import (
        render_rgb_map, make_colorkey_image, sector_corner_pixels,
    )

    sym_key = result.get("sym_key", "D2h")
    ck_n    = 260    # colour-key image side (pixels)

    # ── 1. Pixel-accurate EBSD map ─────────────────────────────────────
    img_rgba, x_min, x_max, y_min, y_max = render_rgb_map(
        result["x"], result["y"],
        result["r"], result["g"], result["b"],
    )
    map_h, map_w = img_rgba.shape[:2]
    dx_map = (x_max - x_min) / max(map_w - 1, 1)
    dy_map = (y_max - y_min) / max(map_h - 1, 1)

    # ── 2. Colour-key triangle image ───────────────────────────────────
    ck_img   = make_colorkey_image(sym_key, color_key=color_key, n=ck_n)
    corners  = sector_corner_pixels(sym_key, n=ck_n)   # (label, col, row, xa, ya)

    # ── 3. Build subplot layout  [map | colour key] ────────────────────
    ck_title = (
        f"IPF-{direction_str} · {sym_label}<br>"
        f"<sup>{'TSL (HKL)' if color_key == 'tsl' else 'HSV (MTEX)'}</sup>"
    )
    fig = make_subplots(
        rows=1, cols=2,
        column_widths=[0.76, 0.24],
        horizontal_spacing=0.04,
        subplot_titles=["", ck_title],
    )

    # ── 4. IPF map trace ───────────────────────────────────────────────
    fig.add_trace(
        go.Image(
            z=img_rgba,
            x0=x_min, dx=dx_map,
            y0=y_min, dy=dy_map,
            hovertemplate="X=%{x:.2f} µm<br>Y=%{y:.2f} µm<extra></extra>",
            name=f"Phase {phase_idx} IPF-{direction_str}",
        ),
        row=1, col=1,
    )
    fig.update_xaxes(
        title_text="X (µm)", range=[x_min, x_max],
        row=1, col=1,
    )
    # Reversed range (y_max first) → y increases downward, matching scan raster
    fig.update_yaxes(
        title_text="Y (µm)",
        range=[y_max, y_min],          # explicit reverse: y_min at top
        scaleanchor="x", scaleratio=1,
        row=1, col=1,
    )

    # ── 5. Colour-key triangle trace ───────────────────────────────────
    _PAD = 8
    fig.add_trace(
        go.Image(
            z=ck_img,
            hoverinfo="skip",
            name="Color key",
        ),
        row=1, col=2,
    )
    fig.update_xaxes(
        showticklabels=False, showgrid=False, zeroline=False,
        range=[-_PAD, ck_n + _PAD],
        row=1, col=2,
    )
    # Reversed range → row 0 at top (standard image convention)
    fig.update_yaxes(
        showticklabels=False, showgrid=False, zeroline=False,
        range=[ck_n + _PAD, -_PAD],   # explicit reverse: row 0 at top
        scaleanchor="x2", scaleratio=1,
        row=1, col=2,
    )

    # ── 6. Corner-direction annotations ───────────────────────────────
    _FONT_CK = dict(size=12, color="#111111", family="Arial")
    for label, col_px, row_px, xanchor, yanchor in corners:
        # Small offset so the text doesn't sit on the coloured pixel
        pad_x = {"left": 5, "right": -5, "center": 0}[xanchor]
        pad_y = {"top": 5, "bottom": -5, "middle": 0}[yanchor]
        fig.add_annotation(
            x=col_px + pad_x,
            y=row_px + pad_y,
            text=f"<b>{label}</b>",
            xref="x2", yref="y2",
            showarrow=False,
            font=_FONT_CK,
            xanchor=xanchor,
            yanchor=yanchor,
            bgcolor="rgba(255,255,255,0.80)",
            bordercolor="rgba(0,0,0,0.25)",
            borderwidth=1,
            borderpad=3,
        )

    # ── 7. Figure-level layout ─────────────────────────────────────────
    fig.update_layout(
        title_text=(
            f"IPF map — Phase {phase_idx}  [{direction_str}]  {sym_label}"
        ),
        title_font_size=14,
        showlegend=False,
        # Transparent background on colour-key subplot so the RGBA alpha
        # channel (outside-sector pixels) shows as plot background
        plot_bgcolor="white",
        paper_bgcolor="#f8f8f8",
        margin=dict(l=60, r=20, t=70, b=50),
    )

    return fig


# ---------------------------------------------------------------------------
# Specimen Reference Frame dialog  (mirrors the MTEX dialog)
# ---------------------------------------------------------------------------

class SpecimenReferenceFrameDialog(QDialog):
    """Modal dialog to choose how to align the CTF spatial frame with
    the Euler (crystal) reference frame — mirroring the MTEX dialog that
    appears when loading EBSD data with ``startupMTEX``.

    The five options are:

    1. **No correction** — use the file as-is (default for data that has
       already been corrected, or when the frame is not important).
    2. **convertEuler2SpatialReferenceFrame** *(MTEX default for CTF/CRC)* —
       adds 180° to every Euler₃ (φ₂) value, rotating each orientation by
       Rz(π) on the specimen side.  This re-expresses orientations from the
       Euler-convention frame (Y up) into the CTF spatial frame (Y down,
       raster-scan convention).
    3. **convertSpatial2EulerReferenceFrame** — flips the Y spatial
       coordinate (y ← y_max − y) so the scan map matches the Euler Y-up
       convention without touching the orientations.
    4. **Custom rotation → Euler angles** — applies an arbitrary Bunge ZXZ
       rotation (entered as φ₁, Φ, φ₂ in degrees) to every pixel orientation
       as a post-multiplication:  R_new = R_pixel · R_custom.
    5. **Custom in-plane rotation → spatial coordinates** — rotates the (X, Y)
       scan map coordinates by an angle (degrees, CCW) around the map centre.
    """

    # Maps internal method key → human-readable label
    _OPTIONS = [
        ("none",            "No correction  (use file as-is)"),
        ("euler2spatial",   "convertEuler2SpatialReferenceFrame  "
                            "(+180° to φ₂) — MTEX default for CTF"),
        ("spatial2euler",   "convertSpatial2EulerReferenceFrame  "
                            "(flip Y coordinate)"),
        ("custom_euler",    "Custom rotation applied to Euler angles"),
        ("custom_spatial",  "Custom in-plane rotation applied to spatial coordinates"),
    ]

    def __init__(self, filename: str = "", parent=None):
        super().__init__(parent)
        self.setWindowTitle("Select Specimen Reference Frame")
        self.setMinimumWidth(560)
        self.setModal(True)

        vbox = QVBoxLayout(self)
        vbox.setSpacing(8)

        # ── header ──────────────────────────────────────────────────────
        title_lbl = QLabel(
            "<b>How should the Euler reference frame relate to the scan "
            "coordinate frame?</b>"
        )
        title_lbl.setWordWrap(True)
        vbox.addWidget(title_lbl)

        if filename:
            fn_lbl = QLabel(f"<small><i>File: {filename}</i></small>")
            fn_lbl.setWordWrap(True)
            vbox.addWidget(fn_lbl)

        sep = QFrame(); sep.setFrameShape(QFrame.HLine); sep.setFrameShadow(QFrame.Sunken)
        vbox.addWidget(sep)

        # ── radio buttons ────────────────────────────────────────────────
        self._btn_group = QButtonGroup(self)
        self._radio_btns: list[QRadioButton] = []

        for i, (key, label) in enumerate(self._OPTIONS):
            rb = QRadioButton(label)
            self._btn_group.addButton(rb, i)
            self._radio_btns.append(rb)
            vbox.addWidget(rb)

        # Default: "euler2spatial" (index 1) — matches MTEX default
        self._radio_btns[1].setChecked(True)

        # ── custom rotation inputs ───────────────────────────────────────
        sep2 = QFrame(); sep2.setFrameShape(QFrame.HLine); sep2.setFrameShadow(QFrame.Sunken)
        vbox.addWidget(sep2)

        self._custom_box = QGroupBox("Custom rotation  (Bunge ZXZ, degrees)")
        custom_form = QFormLayout(self._custom_box)
        custom_form.setSpacing(4)

        self._phi1_spin = QDoubleSpinBox()
        self._phi1_spin.setRange(0.0, 360.0); self._phi1_spin.setDecimals(2)
        self._phi1_spin.setValue(0.0); self._phi1_spin.setSuffix("°")
        custom_form.addRow("φ₁ (Euler1):", self._phi1_spin)

        self._Phi_spin = QDoubleSpinBox()
        self._Phi_spin.setRange(0.0, 180.0); self._Phi_spin.setDecimals(2)
        self._Phi_spin.setValue(0.0); self._Phi_spin.setSuffix("°")
        custom_form.addRow("Φ  (Euler2):", self._Phi_spin)

        self._phi2_spin = QDoubleSpinBox()
        self._phi2_spin.setRange(0.0, 360.0); self._phi2_spin.setDecimals(2)
        self._phi2_spin.setValue(0.0); self._phi2_spin.setSuffix("°")
        custom_form.addRow("φ₂ (Euler3)  /  in-plane angle for spatial:", self._phi2_spin)

        note = QLabel(
            "<small><i>For 'spatial' option only the φ₁ angle is used "
            "(CCW in-plane rotation of the map).</i></small>"
        )
        note.setWordWrap(True)
        custom_form.addRow(note)

        vbox.addWidget(self._custom_box)

        # ── description label ────────────────────────────────────────────
        self._desc_lbl = QLabel()
        self._desc_lbl.setWordWrap(True)
        self._desc_lbl.setStyleSheet("color: #555; font-size: 10px;")
        vbox.addWidget(self._desc_lbl)

        # ── button box ───────────────────────────────────────────────────
        bbox = QDialogButtonBox(
            QDialogButtonBox.Ok | QDialogButtonBox.Cancel, Qt.Horizontal, self
        )
        bbox.accepted.connect(self.accept)
        bbox.rejected.connect(self.reject)
        vbox.addWidget(bbox)

        # Wire radio buttons → show/hide custom box + update description
        self._btn_group.buttonClicked.connect(self._on_radio_changed)
        self._on_radio_changed()   # initialise state

    # ── slots ──────────────────────────────────────────────────────────

    _DESCRIPTIONS = {
        "none": (
            "The orientations and spatial coordinates are used exactly as "
            "stored in the file.  Choose this if your data was already "
            "corrected, or if the reference frame is not important."
        ),
        "euler2spatial": (
            "Adds 180° to every φ₂ (Euler3) angle.  This re-expresses "
            "crystal orientations from the Bunge convention (Y pointing up) "
            "into the Oxford/HKL CTF spatial convention (Y pointing down "
            "in the raster scan).  This is the MTEX default when loading "
            "CTF files."
        ),
        "spatial2euler": (
            "Flips the Y spatial coordinate: y ← y_max − y.  The scan map "
            "is mirrored vertically so that it matches the Euler Y-up "
            "convention.  The Euler angles themselves are not changed."
        ),
        "custom_euler": (
            "Applies the ZXZ Bunge rotation (φ₁, Φ, φ₂) entered above to "
            "every pixel orientation as a post-multiplication: "
            "R_new = R_pixel · R_custom."
        ),
        "custom_spatial": (
            "Rotates the scan map (X, Y) coordinates around the map centre "
            "by φ₁ degrees (counter-clockwise).  Euler angles are not changed."
        ),
    }

    def _on_radio_changed(self, _btn=None):
        method = self.selected_method()
        is_custom = method.startswith("custom_")
        self._custom_box.setVisible(is_custom)
        self._desc_lbl.setText(self._DESCRIPTIONS.get(method, ""))

    # ── public API ──────────────────────────────────────────────────────

    def selected_method(self) -> str:
        idx = self._btn_group.checkedId()
        if 0 <= idx < len(self._OPTIONS):
            return self._OPTIONS[idx][0]
        return "none"

    def custom_rotation(self) -> tuple[float, float, float]:
        """Return (φ₁, Φ, φ₂) in degrees for the custom rotation inputs."""
        return (
            self._phi1_spin.value(),
            self._Phi_spin.value(),
            self._phi2_spin.value(),
        )


# ---------------------------------------------------------------------------
# Worker for VRH averaging
# ---------------------------------------------------------------------------

class _VRHWorker(QThread):
    finished = pyqtSignal(object, float)
    error    = pyqtSignal(str)

    def __init__(self, eb, mb, selected_phases, pressure, temp, method):
        super().__init__()
        self.eb = eb
        self.mb = mb
        self.selected_phases = selected_phases
        self.pressure = pressure
        self.temp = temp
        self.method = method

    def run(self):
        try:
            cij_list, density_list = [], []
            for _, name in self.selected_phases:
                cij = self.mb.get_voigt_matrix_gpa(name, self.pressure, self.temp)
                rho = self.mb.get_density(name)
                if cij is None or rho is None:
                    self.error.emit(f"Cannot load elastic data for '{name}'.")
                    return
                cij_list.append(cij * 1e9)
                density_list.append(rho)
            phase_indices = [idx for idx, _ in self.selected_phases]
            cij_pa, rho = self.eb.vrh_average(
                cij_list, density_list, phase_indices, method=self.method
            )
            if cij_pa is None:
                self.error.emit("VRH averaging returned no result.")
                return
            self.finished.emit(cij_pa / 1e9, rho)
        except Exception as e:
            self.error.emit(str(e))


# ---------------------------------------------------------------------------
# Helper: PlotOptions group box
# ---------------------------------------------------------------------------

_ANIS_SCALARS = ["vp", "vs1", "vs2", "avs", "vpvs1", "vpvs2"]

_SCALAR_LABELS = {
    "vp":    "Vp  — P-wave velocity",
    "vs1":   "Vs1 — fast S-wave",
    "vs2":   "Vs2 — slow S-wave",
    "avs":   "AVs — S-wave anisotropy (%)",
    "vpvs1": "Vp/Vs1 ratio",
    "vpvs2": "Vp/Vs2 ratio",
}


class _PlotOptions(QGroupBox):
    """Reusable plot-options panel with colormap, style, vmin/vmax, point size, colorbar."""

    def __init__(self, title: str = "Plot options",
                 show_cmap: bool = True,
                 default_cmap: str = "Viridis",
                 default_pt_size: int = 3,
                 show_scalar: bool = False,
                 parent=None):
        super().__init__(title, parent)
        f = QFormLayout(self)
        f.setSpacing(3)

        # ── scalar selector (optional) ─────────────────────────────────────
        if show_scalar:
            self.scalar_combo = QComboBox()
            for key in _ANIS_SCALARS:
                self.scalar_combo.addItem(_SCALAR_LABELS.get(key, key), userData=key)
            f.addRow("Display:", self.scalar_combo)
        else:
            self.scalar_combo = None

        # ── stereonet style ───────────────────────────────────────────────
        self.style_combo = QComboBox()
        self.style_combo.addItems(STEREONET_STYLES)
        self.style_combo.setCurrentIndex(0)  # "Smooth fill" default
        f.addRow("Style:", self.style_combo)

        if show_cmap:
            self.cmap_combo = ColorscaleComboBox(default=default_cmap)
            f.addRow("Colormap:", self.cmap_combo)
        else:
            self.cmap_combo = None

        vmin_row = QHBoxLayout()
        self.vmin_spin = QDoubleSpinBox()
        self.vmin_spin.setRange(-1e9, 1e9)
        self.vmin_spin.setValue(0.0)
        self.vmin_spin.setDecimals(3)
        self.vmin_spin.setToolTip("0 = auto")
        self.vmax_spin = QDoubleSpinBox()
        self.vmax_spin.setRange(-1e9, 1e9)
        self.vmax_spin.setValue(0.0)
        self.vmax_spin.setDecimals(3)
        self.vmax_spin.setToolTip("0 = auto")
        vmin_row.addWidget(self.vmin_spin)
        vmin_row.addWidget(QLabel("—"))
        vmin_row.addWidget(self.vmax_spin)
        f.addRow("vmin / vmax:", vmin_row)

        self.pt_size_spin = QSpinBox()
        self.pt_size_spin.setRange(1, 20)
        self.pt_size_spin.setValue(default_pt_size)
        f.addRow("Point size:", self.pt_size_spin)

        self.colorbar_check = QCheckBox("Show colorbar")
        self.colorbar_check.setChecked(True)
        f.addRow(self.colorbar_check)

    @property
    def scalar(self) -> str:
        """Currently selected anisotropy scalar key (e.g. 'vp', 'avs')."""
        if self.scalar_combo is None:
            return "vp"
        return self.scalar_combo.currentData() or self.scalar_combo.currentText()

    @property
    def style(self) -> str:
        return self.style_combo.currentText()

    @property
    def colorscale(self) -> str:
        if self.cmap_combo is None:
            return "Viridis"
        # ColorscaleComboBox exposes .colorscale; plain QComboBox uses .currentText()
        return getattr(self.cmap_combo, "colorscale", self.cmap_combo.currentText())

    @property
    def vmin(self) -> float | None:
        v = self.vmin_spin.value()
        return None if v == 0.0 else v

    @property
    def vmax(self) -> float | None:
        v = self.vmax_spin.value()
        return None if v == 0.0 else v

    @property
    def pt_size(self) -> int:
        return self.pt_size_spin.value()

    @property
    def show_colorbar(self) -> bool:
        return self.colorbar_check.isChecked()


# ---------------------------------------------------------------------------
# EBSD Tab
# ---------------------------------------------------------------------------

class EBSDTab(QWidget):
    def __init__(self, ebsd_backend, material_backend, anisotropy_backend, parent=None):
        super().__init__(parent)
        self.eb = ebsd_backend
        self.mb = material_backend
        self.ab = anisotropy_backend
        self._worker = None
        self._phase_colors: dict[int, str] = {}
        self._phase_color_btns: dict[int, QPushButton] = {}
        self._vrh_grid: dict | None = None
        self._profile_x0 = self._profile_y0 = 0.0
        self._profile_x1 = self._profile_y1 = 100.0
        self._build_ui()

    # ------------------------------------------------------------------
    def _build_ui(self):
        root = QHBoxLayout(self)
        splitter = QSplitter(Qt.Horizontal)
        root.addWidget(splitter)

        # ── Left panel: tabbed control groups ──────────────────────────
        left_scroll = QScrollArea()
        left_scroll.setWidgetResizable(True)
        left_scroll.setFixedWidth(400)
        left_w = QWidget()
        lv = QVBoxLayout(left_w)
        lv.setContentsMargins(2, 2, 2, 2)
        left_scroll.setWidget(left_w)

        # Dataset summary (always visible at top)
        summary_box = QGroupBox("Dataset summary")
        sl = QVBoxLayout(summary_box)
        self.summary_text = QTextEdit()
        self.summary_text.setReadOnly(True)
        mono = QFont("Courier New", 9)
        mono.setStyleHint(QFont.Monospace)
        self.summary_text.setFont(mono)
        self.summary_text.setMinimumHeight(160)
        sl.addWidget(self.summary_text)
        lv.addWidget(summary_box)

        # Tabbed controls
        ctrl_tabs = QTabWidget()
        ctrl_tabs.setTabPosition(QTabWidget.North)
        lv.addWidget(ctrl_tabs)

        # ── Tab 0: Map ────────────────────────────────────────────────
        map_w = QWidget(); map_v = QVBoxLayout(map_w); map_v.setContentsMargins(2, 2, 2, 2)

        self._phase_color_box = QGroupBox("Phase colours")
        self._phase_color_layout = QVBoxLayout(self._phase_color_box)
        self._phase_color_layout.addWidget(QLabel("Load an EBSD file to configure colours."))
        map_v.addWidget(self._phase_color_box)

        pm_opts = _PlotOptions("Phase map options", show_cmap=False, default_pt_size=3)
        self._pm_pt_size = pm_opts.pt_size_spin
        map_v.addWidget(pm_opts)

        roi_box = QGroupBox("Region of Interest (ROI)")
        roi_form = QFormLayout(roi_box)
        roi_form.addRow(QLabel("<small>Hover over the map, read X/Y, enter range.</small>"))
        self.roi_x_min = QDoubleSpinBox(); self.roi_x_min.setRange(-1e9,1e9); self.roi_x_min.setDecimals(1)
        self.roi_x_max = QDoubleSpinBox(); self.roi_x_max.setRange(-1e9,1e9); self.roi_x_max.setDecimals(1)
        self.roi_y_min = QDoubleSpinBox(); self.roi_y_min.setRange(-1e9,1e9); self.roi_y_min.setDecimals(1)
        self.roi_y_max = QDoubleSpinBox(); self.roi_y_max.setRange(-1e9,1e9); self.roi_y_max.setDecimals(1)
        roi_form.addRow("X min:", self.roi_x_min); roi_form.addRow("X max:", self.roi_x_max)
        roi_form.addRow("Y min:", self.roi_y_min); roi_form.addRow("Y max:", self.roi_y_max)
        _rrb = QHBoxLayout()
        self.clip_btn  = QPushButton("Clip to ROI");  self.clip_btn.clicked.connect(self._apply_roi)
        self.reset_btn = QPushButton("Reset ROI");    self.reset_btn.clicked.connect(self._reset_roi)
        _rrb.addWidget(self.clip_btn); _rrb.addWidget(self.reset_btn)
        roi_form.addRow(_rrb)
        map_v.addWidget(roi_box)
        map_v.addStretch()
        ctrl_tabs.addTab(map_w, "Map")

        # ── Tab 1: Pre-process ────────────────────────────────────────
        pp_w = QWidget(); pp_v = QVBoxLayout(pp_w); pp_v.setContentsMargins(2, 2, 2, 2)

        filt_box = QGroupBox("Filter / Threshold")
        filt_form = QFormLayout(filt_box)
        self.mad_spin = QDoubleSpinBox()
        self.mad_spin.setRange(0.0, 5.0); self.mad_spin.setValue(0.7)
        self.mad_spin.setDecimals(2); self.mad_spin.setSuffix("°")
        filt_form.addRow("MAD threshold:", self.mad_spin)
        filt_form.addRow(QPushButton("Apply MAD filter", clicked=self._apply_filter))
        self.downsample_spin = QSpinBox(); self.downsample_spin.setRange(1,100); self.downsample_spin.setValue(10)
        filt_form.addRow("Downsample factor:", self.downsample_spin)
        pp_v.addWidget(filt_box)

        denoise_box = QGroupBox("Denoising (mean orientation filter)")
        dn_form = QFormLayout(denoise_box)
        self.dn_kernel_spin = QSpinBox(); self.dn_kernel_spin.setRange(1,5); self.dn_kernel_spin.setValue(1)
        self.dn_angle_spin  = QDoubleSpinBox(); self.dn_angle_spin.setRange(0.1,30.0); self.dn_angle_spin.setValue(5.0); self.dn_angle_spin.setSuffix("°")
        dn_form.addRow("Kernel size:", self.dn_kernel_spin)
        dn_form.addRow("Max angle:", self.dn_angle_spin)
        dn_form.addRow(QLabel("<small>Replaces each orientation with the quaternion mean\nof same-phase neighbours within the kernel.</small>"))
        dn_form.addRow(QPushButton("Apply denoising", clicked=self._apply_denoise))
        pp_v.addWidget(denoise_box)

        fill_box = QGroupBox("Fill missing data (phase 0 → indexed)")
        fill_form = QFormLayout(fill_box)
        self.fill_method_combo = QComboBox(); self.fill_method_combo.addItems(["nearest", "mean"])
        fill_form.addRow("Method:", self.fill_method_combo)
        fill_form.addRow(QLabel("<small>'nearest' — copy nearest indexed pixel.\n'mean' — quaternion mean of 5 nearest.</small>"))
        fill_form.addRow(QPushButton("Fill missing data", clicked=self._apply_fill))
        pp_v.addWidget(fill_box)

        regrid_box = QGroupBox("Regrid / Resample")
        rg_form = QFormLayout(regrid_box)
        self.rg_mode_combo = QComboBox(); self.rg_mode_combo.addItems(["Scale factor", "Target step (µm)"])
        self.rg_value_spin  = QDoubleSpinBox(); self.rg_value_spin.setRange(0.01, 100.0); self.rg_value_spin.setValue(1.0); self.rg_value_spin.setDecimals(3)
        rg_form.addRow("Mode:", self.rg_mode_combo)
        rg_form.addRow("Value:", self.rg_value_spin)
        rg_form.addRow(QLabel("<small>Scale factor 0.5 = halve step (upsample);\n2.0 = double step (downsample).</small>"))
        rg_form.addRow(QPushButton("Regrid data", clicked=self._apply_regrid))
        pp_v.addWidget(regrid_box)

        sel_box = QGroupBox("Select / Filter by condition")
        sel_form = QFormLayout(sel_box)
        self.sel_phase_combo  = QComboBox(); self.sel_phase_combo.addItem("All phases", userData=None)
        self.sel_mad_max_spin = QDoubleSpinBox(); self.sel_mad_max_spin.setRange(0,10); self.sel_mad_max_spin.setValue(10); self.sel_mad_max_spin.setDecimals(2); self.sel_mad_max_spin.setSuffix("°")
        self.sel_bc_min_spin  = QDoubleSpinBox(); self.sel_bc_min_spin.setRange(0,255); self.sel_bc_min_spin.setValue(0); self.sel_bc_min_spin.setDecimals(0)
        sel_form.addRow("Phase:", self.sel_phase_combo)
        sel_form.addRow("MAD max:", self.sel_mad_max_spin)
        sel_form.addRow("BC min:", self.sel_bc_min_spin)
        sel_form.addRow(QPushButton("Apply selection (clips data)", clicked=self._apply_selection))
        pp_v.addWidget(sel_box)
        pp_v.addStretch()
        ctrl_tabs.addTab(pp_w, "Pre-process")

        # ── Tab 2: Analysis ──────────────────────────────────────────
        an_w = QWidget(); an_v = QVBoxLayout(an_w); an_v.setContentsMargins(2,2,2,2)

        ipf_box = QGroupBox("IPF map (Inverse Pole Figure coloring)")
        ipf_form = QFormLayout(ipf_box)
        self.ipf_phase_combo = QComboBox()
        self.ipf_phase_combo.currentIndexChanged.connect(lambda _: self._update_ipf_sym_label())
        self.ipf_dir_combo   = QComboBox(); self.ipf_dir_combo.addItems(["ND (Z)", "RD (X)", "TD (Y)"])
        # Auto-detected symmetry label (read-only, updated on file load)
        self.ipf_sym_label = QLabel("—")
        self.ipf_sym_label.setStyleSheet("color: #2a7ae2; font-weight: bold;")
        self.ipf_sym_label.setToolTip("Crystal symmetry is auto-detected from the CTF phases table.")
        # Color key selector
        self.ipf_colorkey_combo = QComboBox()
        self.ipf_colorkey_combo.addItem("TSL / HKL  (blue–green–red)", userData="tsl")
        self.ipf_colorkey_combo.addItem("HSV  (white pole, rainbow)", userData="hsv")
        ipf_form.addRow("Phase:", self.ipf_phase_combo)
        ipf_form.addRow("Direction:", self.ipf_dir_combo)
        ipf_form.addRow("Symmetry (auto):", self.ipf_sym_label)
        ipf_form.addRow("Color key:", self.ipf_colorkey_combo)
        ipf_form.addRow(QPushButton("Draw IPF map", clicked=self._draw_ipf_map))
        an_v.addWidget(ipf_box)

        bcbs_box = QGroupBox("Band Contrast / Band Slope map")
        bcbs_form = QFormLayout(bcbs_box)
        self.bcbs_scalar_combo = QComboBox(); self.bcbs_scalar_combo.addItems(["BC", "BS"])
        self.bcbs_cmap_combo   = ColorscaleComboBox(default="Greys")
        self.bcbs_pt_size_spin = QSpinBox(); self.bcbs_pt_size_spin.setRange(1,10); self.bcbs_pt_size_spin.setValue(2)
        bcbs_form.addRow("Scalar:", self.bcbs_scalar_combo)
        bcbs_form.addRow("Colormap:", self.bcbs_cmap_combo)
        bcbs_form.addRow("Point size:", self.bcbs_pt_size_spin)
        bcbs_form.addRow(QPushButton("Draw BC/BS map", clicked=self._draw_bcbs_map))
        an_v.addWidget(bcbs_box)

        kam_box = QGroupBox("KAM  (Kernel Average Misorientation)")
        kam_form = QFormLayout(kam_box)
        self.kam_kernel_spin = QSpinBox(); self.kam_kernel_spin.setRange(1,5); self.kam_kernel_spin.setValue(1)
        self.kam_angle_spin  = QDoubleSpinBox(); self.kam_angle_spin.setRange(0.1,30); self.kam_angle_spin.setValue(5); self.kam_angle_spin.setSuffix("°")
        self.kam_same_phase_chk = QCheckBox("Same phase only"); self.kam_same_phase_chk.setChecked(True)
        self.kam_cmap_combo  = ColorscaleComboBox(default="Hot")
        self.kam_pt_size_spin = QSpinBox(); self.kam_pt_size_spin.setRange(1,10); self.kam_pt_size_spin.setValue(2)
        kam_form.addRow("Kernel size:", self.kam_kernel_spin)
        kam_form.addRow("Max angle:", self.kam_angle_spin)
        kam_form.addRow(self.kam_same_phase_chk)
        kam_form.addRow("Colormap:", self.kam_cmap_combo)
        kam_form.addRow("Point size:", self.kam_pt_size_spin)
        kam_form.addRow(QPushButton("Compute & draw KAM", clicked=self._draw_kam))
        an_v.addWidget(kam_box)

        m2m_box = QGroupBox("Mis2Mean / GROD")
        m2m_form = QFormLayout(m2m_box)
        self.m2m_phase_combo = QComboBox()
        self.m2m_cmap_combo  = ColorscaleComboBox(default="RdBu_r")
        self.m2m_pt_size_spin = QSpinBox(); self.m2m_pt_size_spin.setRange(1,10); self.m2m_pt_size_spin.setValue(2)
        m2m_form.addRow("Phase:", self.m2m_phase_combo)
        m2m_form.addRow("Colormap:", self.m2m_cmap_combo)
        m2m_form.addRow("Point size:", self.m2m_pt_size_spin)
        m2m_form.addRow(QPushButton("Compute & draw Mis2Mean", clicked=self._draw_mis2mean))
        an_v.addWidget(m2m_box)

        euler_box = QGroupBox("Euler-space orientation scatter")
        euler_form = QFormLayout(euler_box)
        self.euler_phase_combo = QComboBox()
        self.euler_ax1_combo   = QComboBox(); self.euler_ax1_combo.addItems(["Euler1","Euler2","Euler3"])
        self.euler_ax2_combo   = QComboBox(); self.euler_ax2_combo.addItems(["Euler2","Euler1","Euler3"]); self.euler_ax2_combo.setCurrentIndex(0)
        self.euler_npts_spin   = QSpinBox(); self.euler_npts_spin.setRange(100,20000); self.euler_npts_spin.setValue(3000); self.euler_npts_spin.setSingleStep(500)
        euler_form.addRow("Phase:", self.euler_phase_combo)
        euler_form.addRow("X axis:", self.euler_ax1_combo)
        euler_form.addRow("Y axis:", self.euler_ax2_combo)
        euler_form.addRow("Max points:", self.euler_npts_spin)
        euler_form.addRow(QPushButton("Draw orientation scatter", clicked=self._draw_euler_scatter))
        an_v.addWidget(euler_box)

        prof_box = QGroupBox("Line profile")
        prof_form = QFormLayout(prof_box)
        self.prof_x0 = QDoubleSpinBox(); self.prof_x0.setRange(-1e9,1e9); self.prof_x0.setDecimals(1)
        self.prof_y0 = QDoubleSpinBox(); self.prof_y0.setRange(-1e9,1e9); self.prof_y0.setDecimals(1)
        self.prof_x1 = QDoubleSpinBox(); self.prof_x1.setRange(-1e9,1e9); self.prof_x1.setDecimals(1)
        self.prof_y1 = QDoubleSpinBox(); self.prof_y1.setRange(-1e9,1e9); self.prof_y1.setDecimals(1)
        self.prof_n_spin = QSpinBox(); self.prof_n_spin.setRange(10,1000); self.prof_n_spin.setValue(100)
        self.prof_scalar_combo = QComboBox()
        self.prof_scalar_combo.addItems(["Euler1","Euler2","Euler3","MAD","BC","BS"])
        prof_form.addRow("Start X:", self.prof_x0); prof_form.addRow("Start Y:", self.prof_y0)
        prof_form.addRow("End X:",   self.prof_x1); prof_form.addRow("End Y:",   self.prof_y1)
        prof_form.addRow("Points:", self.prof_n_spin)
        prof_form.addRow("Scalar:", self.prof_scalar_combo)
        prof_form.addRow(QPushButton("Extract & draw profile", clicked=self._draw_profile))
        an_v.addWidget(prof_box)
        an_v.addStretch()
        ctrl_tabs.addTab(an_w, "Analysis")

        # ── Tab 3: Simulation ─────────────────────────────────────────
        sim_w = QWidget(); sim_v = QVBoxLayout(sim_w); sim_v.setContentsMargins(2,2,2,2)

        sim_box = QGroupBox("Simulate synthetic EBSD")
        sim_form = QFormLayout(sim_box)
        self.sim_ncols_spin = QSpinBox(); self.sim_ncols_spin.setRange(5,500); self.sim_ncols_spin.setValue(50)
        self.sim_nrows_spin = QSpinBox(); self.sim_nrows_spin.setRange(5,500); self.sim_nrows_spin.setValue(50)
        self.sim_step_spin  = QDoubleSpinBox(); self.sim_step_spin.setRange(0.01,100); self.sim_step_spin.setValue(1.0); self.sim_step_spin.setDecimals(3); self.sim_step_spin.setSuffix(" µm")
        self.sim_noise_spin = QDoubleSpinBox(); self.sim_noise_spin.setRange(0,30); self.sim_noise_spin.setValue(0.5); self.sim_noise_spin.setDecimals(2); self.sim_noise_spin.setSuffix("°")
        self.sim_src_combo  = QComboBox(); self.sim_src_combo.addItems(["Random (uniform)", "From current phase 1"])
        sim_form.addRow("Columns:", self.sim_ncols_spin)
        sim_form.addRow("Rows:", self.sim_nrows_spin)
        sim_form.addRow("Step size:", self.sim_step_spin)
        sim_form.addRow("Noise:", self.sim_noise_spin)
        sim_form.addRow("Source:", self.sim_src_combo)
        sim_form.addRow(QPushButton("Generate & load simulation", clicked=self._run_simulation))
        sim_v.addWidget(sim_box)
        sim_v.addStretch()
        ctrl_tabs.addTab(sim_w, "Simulation")

        # ── Tab 4: Anisotropy ─────────────────────────────────────────
        anis_w = QWidget(); anis_v = QVBoxLayout(anis_w); anis_v.setContentsMargins(2,2,2,2)

        self._sn_opts = _PlotOptions("Stereonet options", default_cmap="RdBu_r",
                                     default_pt_size=2, show_scalar=True)
        self._sn_opts.scalar_combo.currentIndexChanged.connect(lambda _: self._replot_stereonet())
        self._sn_opts.style_combo.currentTextChanged.connect(lambda _: self._replot_stereonet())
        self._sn_opts.cmap_combo.currentTextChanged.connect(lambda _: self._replot_stereonet())
        anis_v.addWidget(self._sn_opts)

        self._replot_btn = QPushButton("Replot stereonet")
        self._replot_btn.setToolTip("Redraw with current options — no VRH recalculation")
        self._replot_btn.setEnabled(False)
        self._replot_btn.clicked.connect(self._replot_stereonet)
        anis_v.addWidget(self._replot_btn)

        vrh_box = QGroupBox("Texture averaging (VRH)")
        vrh_form = QFormLayout(vrh_box)
        self.pressure_spin = QDoubleSpinBox(); self.pressure_spin.setRange(0,300); self.pressure_spin.setSuffix(" GPa"); self.pressure_spin.setDecimals(2)
        self.temp_spin     = QDoubleSpinBox(); self.temp_spin.setRange(0,3000); self.temp_spin.setSuffix(" K"); self.temp_spin.setValue(300); self.temp_spin.setDecimals(1)
        self.method_combo  = QComboBox(); self.method_combo.addItems(["voigt","reuss","hill"])
        self.phase_material_table = QTableWidget(); self.phase_material_table.setColumnCount(2)
        self.phase_material_table.setHorizontalHeaderLabels(["EBSD phase", "Material"])
        self.phase_material_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.phase_material_table.setMaximumHeight(150)
        self.vrh_btn = QPushButton("Compute VRH average"); self.vrh_btn.clicked.connect(self._start_vrh)
        vrh_form.addRow("Pressure:", self.pressure_spin)
        vrh_form.addRow("Temperature:", self.temp_spin)
        vrh_form.addRow("Method:", self.method_combo)
        vrh_form.addRow(self.phase_material_table)
        vrh_form.addRow(self.vrh_btn)
        anis_v.addWidget(vrh_box)
        anis_v.addStretch()
        ctrl_tabs.addTab(anis_w, "Anisotropy")

        # ── Tab 5: Export ─────────────────────────────────────────────
        exp_w = QWidget(); exp_v = QVBoxLayout(exp_w); exp_v.setContentsMargins(2,2,2,2)

        exp_box = QGroupBox("Export EBSD data")
        exp_form = QFormLayout(exp_box)
        self.exp_fmt_combo = QComboBox(); self.exp_fmt_combo.addItems(["CSV (.csv)", "CTF (.ctf)", "ANG (.ang)", "HDF5 (.h5)"])
        exp_form.addRow("Format:", self.exp_fmt_combo)
        exp_form.addRow(QPushButton("Export…", clicked=self._export_data))
        exp_v.addWidget(exp_box)
        exp_v.addStretch()
        ctrl_tabs.addTab(exp_w, "Export")

        # Progress bar (shared)
        self.progress = QProgressBar(); self.progress.setRange(0,0); self.progress.setVisible(False)
        lv.addWidget(self.progress)

        splitter.addWidget(left_scroll)

        # ── Right: visualisation tabs ──────────────────────────────────
        right = QWidget()
        rv = QVBoxLayout(right)
        self.vis_tabs = QTabWidget()

        # Phase map
        self.phase_map_plt = PlotlyWidget()
        self.vis_tabs.addTab(self.phase_map_plt, "Phase map")

        # IPF map
        self.ipf_plt = PlotlyWidget()
        self.vis_tabs.addTab(self.ipf_plt, "IPF map")

        # BC/BS map
        self.bcbs_plt = PlotlyWidget()
        self.vis_tabs.addTab(self.bcbs_plt, "BC/BS map")

        # KAM map
        self.kam_plt = PlotlyWidget()
        self.vis_tabs.addTab(self.kam_plt, "KAM")

        # Mis2Mean map
        self.m2m_plt = PlotlyWidget()
        self.vis_tabs.addTab(self.m2m_plt, "Mis2Mean")

        # Orientation scatter (Euler space)
        self.euler_plt = PlotlyWidget()
        self.vis_tabs.addTab(self.euler_plt, "Euler scatter")

        # Line profile
        self.profile_plt = PlotlyWidget()
        self.vis_tabs.addTab(self.profile_plt, "Profile")

        # Pole figure
        self.pf_plt = PlotlyWidget()
        self.vis_tabs.addTab(self.pf_plt, "Pole figure")

        # Anisotropy stereonet
        sn_container = QWidget()
        sn_vbox = QVBoxLayout(sn_container); sn_vbox.setContentsMargins(0,0,0,0)
        sn_btn_row = QHBoxLayout(); sn_btn_row.addStretch()
        self._sn_browser_btn = QPushButton("🌐 Open in browser")
        self._sn_browser_btn.clicked.connect(lambda: self.stereonet_plt.open_in_browser())
        sn_btn_row.addWidget(self._sn_browser_btn)
        sn_vbox.addLayout(sn_btn_row)
        self.stereonet_plt = PlotlyWidget()
        sn_vbox.addWidget(self.stereonet_plt)
        self.vis_tabs.addTab(sn_container, "Anisotropy stereonet")

        # 3-D surface
        self.pv3d = PyVistaWidget()
        self.vis_tabs.addTab(self.pv3d, "3-D surface")

        rv.addWidget(self.vis_tabs)
        splitter.addWidget(right)
        splitter.setSizes([400, 700])

    # ------------------------------------------------------------------
    # Public slots
    # ------------------------------------------------------------------

    def on_file_loaded(self):
        """Called by the main window after a new EBSD file has been loaded.

        Shows the MTEX-style Specimen Reference Frame dialog first so the
        user can choose the appropriate coordinate correction before the
        phase map is drawn or any VRH calculation is attempted.
        """
        self._ask_reference_frame()
        self._update_summary()
        self._refresh_phase_list()
        self._rebuild_phase_color_ui()
        self._refresh_phase_combos()
        self._draw_phase_map()
        self._init_roi_spinboxes()

    def _ask_reference_frame(self):
        """Show the Specimen Reference Frame dialog and apply the correction."""
        filename = getattr(self.eb, 'filename', '') or ''
        import os
        dlg = SpecimenReferenceFrameDialog(
            filename=os.path.basename(filename), parent=self
        )
        if dlg.exec_() == QDialog.Accepted:
            method = dlg.selected_method()
            euler_rot = dlg.custom_rotation() if method.startswith('custom_') else None
            self.eb.apply_reference_frame_correction(method, euler_rot)
            # Invalidate cached VRH grid (orientations may have changed)
            self._vrh_grid = None
            self._replot_btn.setEnabled(False)

    def update_grains_summary(self, grains_backend):
        ebsd_txt  = self.eb.mtex_ebsd_summary()
        grain_txt = self.eb.mtex_grains_summary(grains_backend)
        self.summary_text.setPlainText(ebsd_txt + "\n\n" + grain_txt)

    # ------------------------------------------------------------------
    # Phase colour management
    # ------------------------------------------------------------------

    def _default_color(self, phase_id: int) -> str:
        return DEFAULT_PHASE_COLORS[phase_id % len(DEFAULT_PHASE_COLORS)]

    def _rebuild_phase_color_ui(self):
        """Clear and repopulate the phase colour picker panel."""
        # Remove all existing widgets from layout
        while self._phase_color_layout.count():
            item = self._phase_color_layout.takeAt(0)
            w = item.widget()
            if w:
                w.deleteLater()

        self._phase_color_btns.clear()
        rows = self.eb.phase_rows()   # [(idx, name, pct), ...]

        if not rows:
            self._phase_color_layout.addWidget(QLabel("No phases found."))
            return

        for idx, name, pct in rows:
            if idx not in self._phase_colors:
                self._phase_colors[idx] = self._default_color(idx)
            row_w = QWidget()
            row_l = QHBoxLayout(row_w)
            row_l.setContentsMargins(0, 0, 0, 0)

            lbl = QLabel(f"[{idx}] {name} ({pct:.1f}%)")
            lbl.setMinimumWidth(160)
            row_l.addWidget(lbl)

            btn = QPushButton()
            btn.setFixedSize(44, 22)
            btn.setToolTip("Click to change colour")
            self._set_btn_color(btn, self._phase_colors[idx])
            btn.clicked.connect(lambda _checked, i=idx, b=btn: self._pick_color(i, b))
            row_l.addWidget(btn)
            row_l.addStretch()

            self._phase_color_btns[idx] = btn
            self._phase_color_layout.addWidget(row_w)

        redraw_btn = QPushButton("Redraw phase map")
        redraw_btn.clicked.connect(self._draw_phase_map)
        self._phase_color_layout.addWidget(redraw_btn)

    def _set_btn_color(self, btn: QPushButton, hex_color: str):
        btn.setStyleSheet(
            f"background-color:{hex_color}; border:1px solid #555;"
        )

    def _pick_color(self, phase_id: int, btn: QPushButton):
        current = QColor(self._phase_colors.get(phase_id, "#1f77b4"))
        color = QColorDialog.getColor(current, self, f"Colour for phase {phase_id}")
        if color.isValid():
            self._phase_colors[phase_id] = color.name()
            self._set_btn_color(btn, color.name())
            self._draw_phase_map()

    # ------------------------------------------------------------------
    # Phase map
    # ------------------------------------------------------------------

    def _draw_phase_map(self):
        data = self.eb.map_data()
        if data is None:
            return

        ps   = self._pm_pt_size.value()
        rows = self.eb.phase_rows()   # [(idx, name, pct), ...]
        phase_ids = sorted(set(rows_[0] for rows_ in rows))

        fig = go.Figure()
        for pid, pname, pct in rows:
            mask = data["Phase"].values == pid
            color = self._phase_colors.get(pid, self._default_color(pid))
            fig.add_trace(go.Scattergl(
                x=data["X"].values[mask],
                y=data["Y"].values[mask],
                mode="markers",
                marker=dict(color=color, size=ps, opacity=0.9),
                name=f"[{pid}] {pname} ({pct:.1f}%)",
            ))

        # Mark unindexed (phase == 0)
        mask0 = data["Phase"].values == 0
        if mask0.any():
            fig.add_trace(go.Scattergl(
                x=data["X"].values[mask0],
                y=data["Y"].values[mask0],
                mode="markers",
                marker=dict(color="#cccccc", size=ps, opacity=0.5),
                name="Unindexed",
            ))

        fig.update_layout(
            title="Phase map",
            xaxis_title="X (µm)",
            yaxis_title="Y (µm)",
            yaxis=dict(scaleanchor="x", scaleratio=1),
            legend=dict(itemsizing="constant"),
            hovermode="closest",
        )
        self.phase_map_plt.show_figure(fig)
        self.vis_tabs.setCurrentWidget(self.phase_map_plt)

    # ------------------------------------------------------------------
    # Summary / ROI / filter
    # ------------------------------------------------------------------

    def _update_summary(self):
        self.summary_text.setPlainText(self.eb.mtex_ebsd_summary())

    def _init_roi_spinboxes(self):
        d = self.eb.data
        if d is None:
            return
        self.roi_x_min.setValue(float(d["X"].min()))
        self.roi_x_max.setValue(float(d["X"].max()))
        self.roi_y_min.setValue(float(d["Y"].min()))
        self.roi_y_max.setValue(float(d["Y"].max()))

    def _refresh_phase_list(self):
        rows = self.eb.phase_rows()
        self.phase_material_table.setRowCount(len(rows))
        material_names = self.mb.get_phase_names()
        for row, (idx, name, pct) in enumerate(rows):
            self.phase_material_table.setItem(row, 0, QTableWidgetItem(f"[{idx}] {name}"))
            combo = QComboBox()
            combo.addItems(material_names)
            best = min(range(len(material_names)),
                       key=lambda i, n=name: 0 if material_names[i].lower() == n.lower() else 1,
                       default=0)
            combo.setCurrentIndex(best)
            self.phase_material_table.setCellWidget(row, 1, combo)

    def _apply_filter(self):
        if not self.eb.is_loaded:
            return
        self.eb.filter_mad(self.mad_spin.value())
        self._update_summary()
        self._draw_phase_map()

    def _apply_roi(self):
        if not self.eb.is_loaded:
            return
        self.eb.clip_roi(
            self.roi_x_min.value(), self.roi_x_max.value(),
            self.roi_y_min.value(), self.roi_y_max.value(),
        )
        self._update_summary()
        self._draw_phase_map()

    def _reset_roi(self):
        if not self.eb.is_loaded:
            return
        self.eb.reset_roi()
        self._init_roi_spinboxes()
        self._update_summary()
        self._draw_phase_map()

    # ------------------------------------------------------------------
    # VRH + stereonet
    # ------------------------------------------------------------------

    def _start_vrh(self):
        if not self.eb.is_loaded:
            QMessageBox.warning(self, "No data", "Load an EBSD file first.")
            return
        rows = self.eb.phase_rows()
        selected_phases = []
        for row, (idx, name, _) in enumerate(rows):
            if idx == 0:       # unindexed pixels — must never enter VRH
                continue
            combo = self.phase_material_table.cellWidget(row, 1)
            if combo:
                selected_phases.append((idx, combo.currentText()))
        if not selected_phases:
            return
        self.vrh_btn.setEnabled(False)
        self.progress.setVisible(True)
        self._worker = _VRHWorker(
            self.eb, self.mb, selected_phases,
            self.pressure_spin.value(), self.temp_spin.value(),
            self.method_combo.currentText()
        )
        self._worker.finished.connect(self._on_vrh_done)
        self._worker.error.connect(self._on_vrh_error)
        self._worker.start()

    def _on_vrh_done(self, cij_gpa: np.ndarray, rho: float):
        self.progress.setVisible(False)
        self.vrh_btn.setEnabled(True)
        self.ab.set_from_voigt(cij_gpa * 1e9, rho)

        # Cache the full Cartesian grid — all six scalars are pre-computed inside.
        # _replot_stereonet() reads whichever scalar the user has selected.
        self._vrh_grid = self.ab.compute_stereonet_grid(grid_size=300)
        if self._vrh_grid is not None:
            self._replot_btn.setEnabled(True)
            self._replot_stereonet()

        # 3-D surface (uses whichever scalar the combo shows)
        surf = self.ab.compute_velocity_surface()
        if surf is not None:
            self.pv3d.plot_velocity_surface(surf, scalar=self._sn_opts.scalar)

    def _on_vrh_error(self, msg: str):
        self.progress.setVisible(False)
        self.vrh_btn.setEnabled(True)
        QMessageBox.critical(self, "VRH error", msg)

    def _replot_stereonet(self):
        """Re-render the stereonet from the cached VRH grid — no recalculation."""
        if self._vrh_grid is None:
            return
        opts   = self._sn_opts
        scalar = opts.scalar

        # "Scatter (dots)" needs the raw scattered data, not the Cartesian grid.
        if opts.style == "Scatter (dots)":
            plot_data = self.ab.compute_stereonet_data(step_deg=2.0)
            if plot_data is None:
                return
        else:
            plot_data = self._vrh_grid

        label = _SCALAR_LABELS.get(scalar, scalar.upper())
        self._plot_stereonet(
            plot_data, scalar,
            title=f"{label} — VRH textured aggregate",
        )

        # Keep 3-D surface in sync with selected scalar
        surf = self.ab.compute_velocity_surface()
        if surf is not None:
            self.pv3d.plot_velocity_surface(surf, scalar=scalar)

    def _plot_stereonet(self, data: dict, scalar: str = "vp", title: str = "Stereonet"):
        opts = self._sn_opts

        # Ensure we always pass a grid to the heatmap / contour renderers
        if opts.style != "Scatter (dots)" and "xi" not in data:
            grid_data = self.ab.compute_stereonet_grid(grid_size=300)
            if grid_data is None:
                return
            plot_data = grid_data
        else:
            plot_data = data

        fig = make_stereonet_figure(
            data=plot_data,
            scalar=scalar,
            style=opts.style,
            colorscale=opts.colorscale,
            vmin=opts.vmin,
            vmax=opts.vmax,
            show_colorbar=opts.show_colorbar,
            pt_size=opts.pt_size,
            title=title,
        )
        self.stereonet_plt.show_figure(fig)
        # switch to the stereonet container tab
        for i in range(self.vis_tabs.count()):
            if "stereonet" in self.vis_tabs.tabText(i).lower():
                self.vis_tabs.setCurrentIndex(i)
                break

    # ------------------------------------------------------------------
    # Refresh phase-dependent combo boxes after file load
    # ------------------------------------------------------------------

    def _refresh_phase_combos(self):
        """Re-populate all per-phase combo boxes (IPF, Mis2Mean, etc.)."""
        rows = self.eb.phase_rows()   # [(idx, name, pct), ...]
        combos = [self.ipf_phase_combo, self.m2m_phase_combo,
                  self.euler_phase_combo, self.sel_phase_combo]
        for combo in combos:
            combo.blockSignals(True)
            combo.clear()
            if combo is self.sel_phase_combo:
                combo.addItem("All phases", userData=None)
            for idx, name, _ in rows:
                combo.addItem(f"[{idx}] {name}", userData=idx)
            combo.blockSignals(False)

        # Update auto-detected symmetry label for the first indexed phase
        self._update_ipf_sym_label()

    # ------------------------------------------------------------------
    # IPF map
    # ------------------------------------------------------------------

    def _update_ipf_sym_label(self):
        """Refresh the auto-detected symmetry label from the current IPF phase."""
        if not self.eb.is_loaded:
            self.ipf_sym_label.setText("—")
            return
        phase_idx = self.ipf_phase_combo.currentData()
        if phase_idx is None:
            self.ipf_sym_label.setText("—")
            return
        try:
            label = self.eb.phase_symmetry_label(phase_idx)
            self.ipf_sym_label.setText(label)
        except Exception:
            self.ipf_sym_label.setText("—")

    def _draw_ipf_map(self):
        if not self.eb.is_loaded:
            QMessageBox.warning(self, "No data", "Load an EBSD file first.")
            return
        phase_idx = self.ipf_phase_combo.currentData()
        if phase_idx is None:
            QMessageBox.warning(self, "No phase", "No phase selected.")
            return
        direction_str = self.ipf_dir_combo.currentText().split()[0]   # "ND", "RD", "TD"
        color_key = self.ipf_colorkey_combo.currentData() or "tsl"

        # Refresh label in case the phase combo changed without triggering update
        self._update_ipf_sym_label()

        self.progress.setVisible(True)
        try:
            result = self.eb.ipf_map_colors(phase_idx, direction_str, color_key)
        except Exception as e:
            QMessageBox.critical(self, "IPF error", str(e))
            self.progress.setVisible(False)
            return
        self.progress.setVisible(False)

        if result is None:
            QMessageBox.information(self, "IPF map", "No data for selected phase.")
            return

        sym_label = self.ipf_sym_label.text()

        try:
            fig = _build_ipf_figure(result, phase_idx, direction_str,
                                    color_key, sym_label)
        except Exception as e:
            # Fallback to plain scatter
            fig = go.Figure()
            fig.add_trace(go.Scattergl(
                x=result["x"], y=result["y"],
                mode="markers",
                marker=dict(color=result["rgb_hex"], size=3, opacity=0.9),
                hovertemplate="X=%{x:.1f}<br>Y=%{y:.1f}<extra></extra>",
                name=f"Phase {phase_idx} IPF-{direction_str}",
            ))
            fig.update_layout(
                title=f"IPF map — Phase {phase_idx}, direction {direction_str}  "
                      f"[fallback: {e}]",
                xaxis_title="X (µm)", yaxis_title="Y (µm)",
                yaxis=dict(scaleanchor="x", scaleratio=1),
            )

        self.ipf_plt.show_figure(fig)
        self._switch_to_tab("IPF map")

    # ------------------------------------------------------------------
    # BC / BS map
    # ------------------------------------------------------------------

    def _draw_bcbs_map(self):
        if not self.eb.is_loaded:
            QMessageBox.warning(self, "No data", "Load an EBSD file first.")
            return
        result = self.eb.bc_bs_map_data()
        if result is None:
            QMessageBox.information(self, "BC/BS", "No BC/BS data in file.")
            return
        scalar_key = self.bcbs_scalar_combo.currentText().lower()
        cmap = self.bcbs_cmap_combo.currentText()
        ps   = self.bcbs_pt_size_spin.value()
        vals = result.get(scalar_key)
        if vals is None:
            QMessageBox.information(self, "BC/BS", f"Column '{scalar_key.upper()}' not found.")
            return
        fig = go.Figure()
        fig.add_trace(go.Scattergl(
            x=result["x"], y=result["y"],
            mode="markers",
            marker=dict(color=vals, colorscale=cmap, showscale=True,
                        colorbar=dict(title=scalar_key.upper()), size=ps, opacity=0.9),
            hovertemplate=f"X=%{{x:.1f}}<br>Y=%{{y:.1f}}<br>{scalar_key.upper()}=%{{marker.color:.0f}}<extra></extra>",
        ))
        fig.update_layout(
            title=f"{scalar_key.upper()} map",
            xaxis_title="X (µm)", yaxis_title="Y (µm)",
            yaxis=dict(scaleanchor="x", scaleratio=1),
        )
        self.bcbs_plt.show_figure(fig)
        self._switch_to_tab("BC/BS map")

    # ------------------------------------------------------------------
    # KAM
    # ------------------------------------------------------------------

    def _draw_kam(self):
        if not self.eb.is_loaded:
            QMessageBox.warning(self, "No data", "Load an EBSD file first.")
            return
        self.progress.setVisible(True)
        try:
            kam = self.eb.compute_kam(
                kernel_size=self.kam_kernel_spin.value(),
                max_angle_deg=self.kam_angle_spin.value(),
                same_phase_only=self.kam_same_phase_chk.isChecked(),
            )
        except Exception as e:
            QMessageBox.critical(self, "KAM error", str(e))
            self.progress.setVisible(False)
            return
        self.progress.setVisible(False)
        if kam is None:
            return
        data = self.eb.map_data()
        valid = ~np.isnan(kam)
        cmap = self.kam_cmap_combo.currentText()
        ps   = self.kam_pt_size_spin.value()
        fig = go.Figure()
        fig.add_trace(go.Scattergl(
            x=data["X"].values[valid], y=data["Y"].values[valid],
            mode="markers",
            marker=dict(color=kam[valid], colorscale=cmap, showscale=True,
                        colorbar=dict(title="KAM (°)"), size=ps, opacity=0.9),
            hovertemplate="X=%{x:.1f}<br>Y=%{y:.1f}<br>KAM=%{marker.color:.3f}°<extra></extra>",
        ))
        fig.update_layout(
            title="KAM — Kernel Average Misorientation",
            xaxis_title="X (µm)", yaxis_title="Y (µm)",
            yaxis=dict(scaleanchor="x", scaleratio=1),
        )
        self.kam_plt.show_figure(fig)
        self._switch_to_tab("KAM")

    # ------------------------------------------------------------------
    # Mis2Mean
    # ------------------------------------------------------------------

    def _draw_mis2mean(self):
        if not self.eb.is_loaded:
            QMessageBox.warning(self, "No data", "Load an EBSD file first.")
            return
        phase_idx = self.m2m_phase_combo.currentData()
        if phase_idx is None:
            return
        self.progress.setVisible(True)
        try:
            m2m = self.eb.compute_mis2mean(phase_idx)
        except Exception as e:
            QMessageBox.critical(self, "Mis2Mean error", str(e))
            self.progress.setVisible(False)
            return
        self.progress.setVisible(False)
        if m2m is None:
            return
        data = self.eb.map_data()
        valid = ~np.isnan(m2m)
        cmap = self.m2m_cmap_combo.currentText()
        ps   = self.m2m_pt_size_spin.value()
        fig = go.Figure()
        fig.add_trace(go.Scattergl(
            x=data["X"].values[valid], y=data["Y"].values[valid],
            mode="markers",
            marker=dict(color=m2m[valid], colorscale=cmap, showscale=True,
                        colorbar=dict(title="Mis2Mean (°)"), size=ps, opacity=0.9),
            hovertemplate="X=%{x:.1f}<br>Y=%{y:.1f}<br>Mis2Mean=%{marker.color:.3f}°<extra></extra>",
        ))
        fig.update_layout(
            title=f"Mis2Mean — Phase {phase_idx}",
            xaxis_title="X (µm)", yaxis_title="Y (µm)",
            yaxis=dict(scaleanchor="x", scaleratio=1),
        )
        self.m2m_plt.show_figure(fig)
        self._switch_to_tab("Mis2Mean")

    # ------------------------------------------------------------------
    # Euler-space orientation scatter
    # ------------------------------------------------------------------

    def _draw_euler_scatter(self):
        if not self.eb.is_loaded:
            QMessageBox.warning(self, "No data", "Load an EBSD file first.")
            return
        phase_idx  = self.euler_phase_combo.currentData()
        ax1        = self.euler_ax1_combo.currentText()
        ax2        = self.euler_ax2_combo.currentText()
        max_pts    = self.euler_npts_spin.value()
        result = self.eb.orientation_scatter_data(phase_idx, (ax1, ax2), max_pts)
        if result is None:
            return
        fig = go.Figure()
        fig.add_trace(go.Scattergl(
            x=result.get(ax1, []), y=result.get(ax2, []),
            mode="markers",
            marker=dict(size=3, color=result.get("Euler2", result.get(ax1)),
                        colorscale="Viridis", showscale=True,
                        colorbar=dict(title="Φ (°)")),
            hovertemplate=f"{ax1}=%{{x:.1f}}°<br>{ax2}=%{{y:.1f}}°<extra></extra>",
        ))
        fig.update_layout(
            title=f"Orientation scatter — Phase {phase_idx}",
            xaxis_title=f"{ax1} (°)", yaxis_title=f"{ax2} (°)",
        )
        self.euler_plt.show_figure(fig)
        self._switch_to_tab("Euler scatter")

    # ------------------------------------------------------------------
    # Line profile
    # ------------------------------------------------------------------

    def _draw_profile(self):
        if not self.eb.is_loaded:
            QMessageBox.warning(self, "No data", "Load an EBSD file first.")
            return
        x0 = self.prof_x0.value(); y0 = self.prof_y0.value()
        x1 = self.prof_x1.value(); y1 = self.prof_y1.value()
        n  = self.prof_n_spin.value()
        col = self.prof_scalar_combo.currentText()

        profile = self.eb.get_line_profile(x0, y0, x1, y1, n, scalars=[col])
        if profile is None or col not in profile:
            QMessageBox.information(self, "Profile", f"Column '{col}' not found.")
            return

        fig = go.Figure()
        # Scalar plot along the profile
        fig.add_trace(go.Scatter(
            x=profile["distance"], y=profile[col],
            mode="lines+markers",
            marker=dict(size=4),
            name=col,
            hovertemplate="d=%{x:.1f} µm<br>" + col + "=%{y:.2f}<extra></extra>",
        ))
        fig.update_layout(
            title=f"Line profile  ({x0:.0f},{y0:.0f}) → ({x1:.0f},{y1:.0f})",
            xaxis_title="Distance (µm)",
            yaxis_title=col,
        )
        self.profile_plt.show_figure(fig)
        self._switch_to_tab("Profile")

        # Overlay the profile line on the phase map
        data = self.eb.map_data()
        if data is not None:
            import plotly.graph_objects as _go
            fig2 = self.phase_map_plt._fig if hasattr(self.phase_map_plt, '_fig') else None
            # re-draw phase map with the line overlaid
            self._draw_phase_map()
            # After redraw, add a shape for the line
            self.phase_map_plt.fig_handle.add_shape(
                type="line", xref="x", yref="y",
                x0=x0, y0=y0, x1=x1, y1=y1,
                line=dict(color="red", width=2, dash="dash"),
            ) if hasattr(self.phase_map_plt, 'fig_handle') else None

    # ------------------------------------------------------------------
    # Pre-processing actions
    # ------------------------------------------------------------------

    def _apply_denoise(self):
        if not self.eb.is_loaded:
            QMessageBox.warning(self, "No data", "Load an EBSD file first.")
            return
        self.progress.setVisible(True)
        try:
            self.eb.denoise(
                kernel_size=self.dn_kernel_spin.value(),
                max_angle_deg=self.dn_angle_spin.value(),
            )
        except Exception as e:
            QMessageBox.critical(self, "Denoise error", str(e))
        finally:
            self.progress.setVisible(False)
        self._update_summary()
        self._draw_phase_map()

    def _apply_fill(self):
        if not self.eb.is_loaded:
            QMessageBox.warning(self, "No data", "Load an EBSD file first.")
            return
        self.progress.setVisible(True)
        try:
            self.eb.fill_missing(self.fill_method_combo.currentText())
        except Exception as e:
            QMessageBox.critical(self, "Fill error", str(e))
        finally:
            self.progress.setVisible(False)
        self._update_summary()
        self._draw_phase_map()

    def _apply_regrid(self):
        if not self.eb.is_loaded:
            QMessageBox.warning(self, "No data", "Load an EBSD file first.")
            return
        self.progress.setVisible(True)
        try:
            val = self.rg_value_spin.value()
            if self.rg_mode_combo.currentIndex() == 0:
                self.eb.regrid_data(scale=val)
            else:
                self.eb.regrid_data(target_step=val)
        except Exception as e:
            QMessageBox.critical(self, "Regrid error", str(e))
        finally:
            self.progress.setVisible(False)
        self._update_summary()
        self._init_roi_spinboxes()
        self._draw_phase_map()

    def _apply_selection(self):
        if not self.eb.is_loaded:
            QMessageBox.warning(self, "No data", "Load an EBSD file first.")
            return
        phase_data = self.sel_phase_combo.currentData()
        phase_val  = int(phase_data) if phase_data is not None else None
        mad_max    = self.sel_mad_max_spin.value()
        bc_min     = self.sel_bc_min_spin.value()

        kwargs: dict = {}
        if phase_val is not None:
            kwargs["phase"] = phase_val
        if mad_max < 10.0:
            kwargs["mad_max"] = mad_max
        if bc_min > 0:
            kwargs["bc_min"] = bc_min

        filtered = self.eb.select_by_condition(**kwargs)
        if filtered is None or len(filtered) == 0:
            QMessageBox.information(self, "Select", "No pixels match the condition.")
            return
        # Replace active data with selection
        self.eb.data = filtered
        self._update_summary()
        self._draw_phase_map()

    # ------------------------------------------------------------------
    # Simulation
    # ------------------------------------------------------------------

    def _run_simulation(self):
        import numpy as _np
        n_cols = self.sim_ncols_spin.value()
        n_rows = self.sim_nrows_spin.value()
        step   = self.sim_step_spin.value()
        noise  = self.sim_noise_spin.value()

        if self.sim_src_combo.currentIndex() == 0:
            # Uniform random orientations
            rng = _np.random.default_rng(42)
            eulers = rng.uniform([0, 0, 0], [360, 90, 360], size=(200, 3))
        else:
            # Use current phase-1 orientations
            if not self.eb.is_loaded:
                QMessageBox.warning(self, "No data", "Load an EBSD file first.")
                return
            ph1 = self.eb.euler_angles(1)
            if ph1 is None or len(ph1) == 0:
                QMessageBox.warning(self, "No data", "No phase-1 orientations found.")
                return
            eulers = ph1[["Euler1", "Euler2", "Euler3"]].to_numpy()

        df = self.eb.simulate_from_euler(eulers, n_cols, n_rows, step, noise_deg=noise)
        self.eb.load_synthetic(df)
        self.on_file_loaded()
        QMessageBox.information(self, "Simulation",
                                f"Generated {n_cols}×{n_rows} = {len(df)} pixel synthetic map.")

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def _export_data(self):
        if not self.eb.is_loaded:
            QMessageBox.warning(self, "No data", "Load an EBSD file first.")
            return
        from PyQt5.QtWidgets import QFileDialog
        fmt_text = self.exp_fmt_combo.currentText()
        fmt_map = {
            "CSV (.csv)": ("csv", "CSV files (*.csv)"),
            "CTF (.ctf)": ("ctf", "CTF files (*.ctf)"),
            "ANG (.ang)": ("ang", "ANG files (*.ang)"),
            "HDF5 (.h5)": ("hdf5", "HDF5 files (*.h5 *.hdf5)"),
        }
        fmt_key, file_filter = fmt_map.get(fmt_text, ("csv", "CSV files (*.csv)"))
        path, _ = QFileDialog.getSaveFileName(self, "Export EBSD data", "", file_filter)
        if not path:
            return
        try:
            self.eb.export(path, fmt=fmt_key)
            QMessageBox.information(self, "Export", f"Saved to:\n{path}")
        except Exception as e:
            QMessageBox.critical(self, "Export error", str(e))

    # ------------------------------------------------------------------
    # Helper: switch visible tab
    # ------------------------------------------------------------------

    def _switch_to_tab(self, tab_name: str):
        """Activate the visualisation tab whose text contains *tab_name*."""
        for i in range(self.vis_tabs.count()):
            if tab_name.lower() in self.vis_tabs.tabText(i).lower():
                self.vis_tabs.setCurrentIndex(i)
                return
