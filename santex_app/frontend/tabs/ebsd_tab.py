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
    QScrollArea, QCheckBox,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QFont, QColor

from frontend.widgets.plotly_widget import PlotlyWidget, COLORSCALES, DEFAULT_PHASE_COLORS
from frontend.tabs._stereonet import make_stereonet_figure, STEREONET_STYLES
from frontend.widgets.pyvista_widget import PyVistaWidget

import plotly.graph_objects as go


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
            self.cmap_combo = QComboBox()
            self.cmap_combo.addItems(COLORSCALES)
            idx = COLORSCALES.index(default_cmap) if default_cmap in COLORSCALES else 0
            self.cmap_combo.setCurrentIndex(idx)
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
        return self.cmap_combo.currentText() if self.cmap_combo else "Viridis"

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
        self._phase_colors: dict[int, str] = {}   # phase_id → hex color
        self._phase_color_btns: dict[int, QPushButton] = {}
        self._vrh_grid: dict | None = None          # cached stereonet grid after VRH
        self._build_ui()

    # ------------------------------------------------------------------
    def _build_ui(self):
        root = QHBoxLayout(self)
        splitter = QSplitter(Qt.Horizontal)
        root.addWidget(splitter)

        # ---- Left: scrollable controls ----
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFixedWidth(390)
        left = QWidget()
        lv = QVBoxLayout(left)
        lv.setContentsMargins(4, 4, 4, 4)
        scroll.setWidget(left)

        # Dataset summary
        summary_box = QGroupBox("Dataset summary")
        sl = QVBoxLayout(summary_box)
        self.summary_text = QTextEdit()
        self.summary_text.setReadOnly(True)
        mono = QFont("Courier New", 9)
        mono.setStyleHint(QFont.Monospace)
        self.summary_text.setFont(mono)
        self.summary_text.setMinimumHeight(180)
        sl.addWidget(self.summary_text)
        lv.addWidget(summary_box)

        # Phase colour pickers
        self._phase_color_box = QGroupBox("Phase colours")
        self._phase_color_layout = QVBoxLayout(self._phase_color_box)
        self._phase_color_layout.addWidget(QLabel("Load an EBSD file to configure colours."))
        lv.addWidget(self._phase_color_box)

        # Phase map options
        pm_opts = _PlotOptions("Phase map options", show_cmap=False, default_pt_size=3)
        self._pm_pt_size = pm_opts.pt_size_spin
        lv.addWidget(pm_opts)

        # ROI
        roi_box = QGroupBox("Region of Interest (ROI)")
        roi_form = QFormLayout(roi_box)
        roi_form.addRow(QLabel(
            "<small>Hover over the map to read X/Y, then enter the range here.</small>"
        ))
        self.roi_x_min = QDoubleSpinBox()
        self.roi_x_min.setRange(-1e9, 1e9); self.roi_x_min.setDecimals(1)
        self.roi_x_max = QDoubleSpinBox()
        self.roi_x_max.setRange(-1e9, 1e9); self.roi_x_max.setDecimals(1)
        self.roi_y_min = QDoubleSpinBox()
        self.roi_y_min.setRange(-1e9, 1e9); self.roi_y_min.setDecimals(1)
        self.roi_y_max = QDoubleSpinBox()
        self.roi_y_max.setRange(-1e9, 1e9); self.roi_y_max.setDecimals(1)
        roi_form.addRow("X min:", self.roi_x_min)
        roi_form.addRow("X max:", self.roi_x_max)
        roi_form.addRow("Y min:", self.roi_y_min)
        roi_form.addRow("Y max:", self.roi_y_max)
        roi_btn_row = QHBoxLayout()
        self.clip_btn  = QPushButton("Clip to ROI")
        self.clip_btn.clicked.connect(self._apply_roi)
        self.reset_btn = QPushButton("Reset ROI")
        self.reset_btn.clicked.connect(self._reset_roi)
        roi_btn_row.addWidget(self.clip_btn)
        roi_btn_row.addWidget(self.reset_btn)
        roi_form.addRow(roi_btn_row)
        lv.addWidget(roi_box)

        # Pre-processing
        filt_box = QGroupBox("Pre-processing")
        filt_form = QFormLayout(filt_box)
        self.mad_spin = QDoubleSpinBox()
        self.mad_spin.setRange(0.0, 5.0)
        self.mad_spin.setValue(0.7)
        self.mad_spin.setDecimals(2)
        self.mad_spin.setSuffix("°")
        filt_form.addRow("MAD threshold:", self.mad_spin)
        filt_form.addRow(QPushButton("Apply MAD filter",
                                     clicked=self._apply_filter))
        self.downsample_spin = QSpinBox()
        self.downsample_spin.setRange(1, 100)
        self.downsample_spin.setValue(10)
        filt_form.addRow("Downsample factor:", self.downsample_spin)
        lv.addWidget(filt_box)

        # Stereonet options — scalar selector lives here so users can switch
        # metrics and replot *without* rerunning the VRH calculation.
        self._sn_opts = _PlotOptions("Stereonet options", default_cmap="RdBu_r",
                                     default_pt_size=2, show_scalar=True)
        # Auto-replot when scalar, style or colormap changes (only if VRH is cached)
        self._sn_opts.scalar_combo.currentIndexChanged.connect(
            lambda _: self._replot_stereonet()
        )
        self._sn_opts.style_combo.currentTextChanged.connect(
            lambda _: self._replot_stereonet()
        )
        self._sn_opts.cmap_combo.currentTextChanged.connect(
            lambda _: self._replot_stereonet()
        )
        lv.addWidget(self._sn_opts)

        self._replot_btn = QPushButton("Replot stereonet")
        self._replot_btn.setToolTip(
            "Redraw the stereonet with the current style / colormap / scalar.\n"
            "No VRH recalculation — uses the last computed result."
        )
        self._replot_btn.setEnabled(False)
        self._replot_btn.clicked.connect(self._replot_stereonet)
        lv.addWidget(self._replot_btn)

        # VRH
        vrh_box = QGroupBox("Texture averaging (VRH)")
        vrh_form = QFormLayout(vrh_box)
        self.pressure_spin = QDoubleSpinBox()
        self.pressure_spin.setRange(0, 300)
        self.pressure_spin.setSuffix(" GPa")
        self.pressure_spin.setDecimals(2)
        vrh_form.addRow("Pressure:", self.pressure_spin)
        self.temp_spin = QDoubleSpinBox()
        self.temp_spin.setRange(0, 3000)
        self.temp_spin.setSuffix(" K")
        self.temp_spin.setValue(300)
        self.temp_spin.setDecimals(1)
        vrh_form.addRow("Temperature:", self.temp_spin)
        self.method_combo = QComboBox()
        self.method_combo.addItems(["voigt", "reuss", "hill"])
        vrh_form.addRow("Method:", self.method_combo)
        self.phase_material_table = QTableWidget()
        self.phase_material_table.setColumnCount(2)
        self.phase_material_table.setHorizontalHeaderLabels(["EBSD phase", "Material"])
        self.phase_material_table.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.phase_material_table.setMaximumHeight(150)
        vrh_form.addRow(self.phase_material_table)
        self.vrh_btn = QPushButton("Compute VRH average")
        self.vrh_btn.clicked.connect(self._start_vrh)
        vrh_form.addRow(self.vrh_btn)
        lv.addWidget(vrh_box)

        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setVisible(False)
        lv.addWidget(self.progress)
        lv.addStretch()

        splitter.addWidget(scroll)

        # ---- Right: visualisation tabs ----
        right = QWidget()
        rv = QVBoxLayout(right)
        self.vis_tabs = QTabWidget()

        self.phase_map_plt = PlotlyWidget()
        self.vis_tabs.addTab(self.phase_map_plt, "Phase map")

        self.pf_plt = PlotlyWidget()
        self.vis_tabs.addTab(self.pf_plt, "Pole figure")

        # Stereonet tab — wrap in a widget so we can add the "Open in browser" button
        sn_container = QWidget()
        sn_vbox = QVBoxLayout(sn_container)
        sn_vbox.setContentsMargins(0, 0, 0, 0)
        sn_btn_row = QHBoxLayout()
        sn_btn_row.addStretch()
        self._sn_browser_btn = QPushButton("🌐 Open in browser")
        self._sn_browser_btn.setToolTip("Open this stereonet in the system browser for a larger interactive view")
        self._sn_browser_btn.clicked.connect(lambda: self.stereonet_plt.open_in_browser())
        sn_btn_row.addWidget(self._sn_browser_btn)
        sn_vbox.addLayout(sn_btn_row)
        self.stereonet_plt = PlotlyWidget()
        sn_vbox.addWidget(self.stereonet_plt)
        self.vis_tabs.addTab(sn_container, "Anisotropy stereonet")

        self.pv3d = PyVistaWidget()
        self.vis_tabs.addTab(self.pv3d, "3-D surface")

        rv.addWidget(self.vis_tabs)
        splitter.addWidget(right)
        splitter.setSizes([390, 710])

    # ------------------------------------------------------------------
    # Public slots
    # ------------------------------------------------------------------

    def on_file_loaded(self):
        self._update_summary()
        self._refresh_phase_list()
        self._rebuild_phase_color_ui()
        self._draw_phase_map()
        self._init_roi_spinboxes()

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
