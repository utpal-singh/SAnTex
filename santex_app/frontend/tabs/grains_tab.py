"""
GrainsTab — grain analysis with Plotly visualisation.

Sub-tabs
--------
  Reconstruction   : grain map, grain size map
  Shape            : histograms (area, aspect ratio, GOS), scatter size vs shape
  Orientation      : KAM map, GROD map
  Boundaries       : boundary map (scatter at midpoints), misorientation dist, CSL, twist/tilt
  Triple Points    : triple-point map
  Export           : CSV / Neper / CTF
"""

from __future__ import annotations
import numpy as np

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QDoubleSpinBox, QSpinBox, QPushButton, QSplitter,
    QComboBox, QFormLayout, QTabWidget, QProgressBar,
    QTableWidget, QTableWidgetItem, QHeaderView, QAbstractItemView,
    QCheckBox, QFileDialog, QMessageBox, QTextEdit, QSizePolicy,
    QScrollArea,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal

from frontend.widgets.plotly_widget import PlotlyWidget, COLORSCALES, DEFAULT_PHASE_COLORS
from frontend.tabs.ebsd_tab import _PlotOptions

import plotly.graph_objects as go
from plotly.subplots import make_subplots


# ---------------------------------------------------------------------------
# Workers
# ---------------------------------------------------------------------------

class _ReconWorker(QThread):
    finished = pyqtSignal()
    error    = pyqtSignal(str)

    def __init__(self, gb, data, threshold, min_pixels):
        super().__init__()
        self.gb = gb; self.data = data
        self.threshold = threshold; self.min_pixels = min_pixels

    def run(self):
        try:
            self.gb.reconstruct(self.data, self.threshold, self.min_pixels)
            self.finished.emit()
        except Exception as e:
            self.error.emit(str(e))


class _AnalysisWorker(QThread):
    finished = pyqtSignal(object)
    error    = pyqtSignal(str)

    def __init__(self, fn):
        super().__init__()
        self.fn = fn

    def run(self):
        try:
            self.finished.emit(self.fn())
        except Exception as e:
            self.error.emit(str(e))


# ---------------------------------------------------------------------------
# Main tab
# ---------------------------------------------------------------------------

class GrainsTab(QWidget):
    def __init__(self, grains_backend, ebsd_backend, ebsd_tab=None, parent=None):
        super().__init__(parent)
        self.gb = grains_backend
        self.eb = ebsd_backend
        self.ebsd_tab = ebsd_tab
        self._worker = None
        self._build_ui()

    def _build_ui(self):
        root = QVBoxLayout(self)
        self.sub_tabs = QTabWidget()
        root.addWidget(self.sub_tabs)
        self.sub_tabs.addTab(self._build_recon_panel(),    "Reconstruction")
        self.sub_tabs.addTab(self._build_shape_panel(),    "Shape Parameters")
        self.sub_tabs.addTab(self._build_orient_panel(),   "Orientation")
        self.sub_tabs.addTab(self._build_boundary_panel(), "Grain Boundaries")
        self.sub_tabs.addTab(self._build_tp_panel(),       "Triple Points")
        self.sub_tabs.addTab(self._build_export_panel(),   "Export")

    # ==================================================================
    # Reconstruction
    # ==================================================================

    def _build_recon_panel(self):
        panel = QWidget()
        layout = QHBoxLayout(panel)
        splitter = QSplitter(Qt.Horizontal)
        layout.addWidget(splitter)

        left = QWidget(); lv = QVBoxLayout(left)
        ctrl = QGroupBox("Reconstruction parameters")
        cf = QFormLayout(ctrl)
        self.threshold_spin = QDoubleSpinBox()
        self.threshold_spin.setRange(0.1, 90.0)
        self.threshold_spin.setSuffix("°")
        self.threshold_spin.setDecimals(1)
        self.threshold_spin.setValue(10.0)
        cf.addRow("GB angle threshold:", self.threshold_spin)
        self.minpx_spin = QSpinBox()
        self.minpx_spin.setRange(1, 10000)
        self.minpx_spin.setValue(2)
        cf.addRow("Min pixels per grain:", self.minpx_spin)
        lv.addWidget(ctrl)

        merge_box = QGroupBox("Post-processing")
        mf = QFormLayout(merge_box)
        self.merge_check = QCheckBox("Merge grains below:")
        self.merge_angle_spin = QDoubleSpinBox()
        self.merge_angle_spin.setRange(0.1, 30.0)
        self.merge_angle_spin.setSuffix("°")
        self.merge_angle_spin.setValue(5.0)
        mf.addRow(self.merge_check, self.merge_angle_spin)
        lv.addWidget(merge_box)

        self._recon_map_opts = _PlotOptions("Map options", default_cmap="Turbo",
                                             default_pt_size=2)
        lv.addWidget(self._recon_map_opts)

        self.run_btn = QPushButton("Run reconstruction")
        self.run_btn.clicked.connect(self._run_recon)
        lv.addWidget(self.run_btn)
        self.recon_progress = QProgressBar()
        self.recon_progress.setRange(0, 0)
        self.recon_progress.setVisible(False)
        lv.addWidget(self.recon_progress)
        self.recon_info = QTextEdit()
        self.recon_info.setReadOnly(True)
        self.recon_info.setMaximumHeight(100)
        lv.addWidget(self.recon_info)
        lv.addStretch()
        splitter.addWidget(left)

        right = QWidget(); rv = QVBoxLayout(right)
        self.recon_vis_tabs = QTabWidget()
        self.grain_map_plt = PlotlyWidget()
        self.size_map_plt  = PlotlyWidget()
        self.recon_vis_tabs.addTab(self.grain_map_plt, "Grain ID map")
        self.recon_vis_tabs.addTab(self.size_map_plt,  "Grain size map")
        rv.addWidget(self.recon_vis_tabs)
        splitter.addWidget(right)
        splitter.setSizes([340, 760])
        return panel

    def _run_recon(self):
        if not self.eb.is_loaded:
            QMessageBox.warning(self, "No EBSD data", "Load an EBSD file first.")
            return
        data = self.eb.data
        if data is None:
            return
        self.run_btn.setEnabled(False)
        self.recon_progress.setVisible(True)
        self._worker = _ReconWorker(
            self.gb, data, self.threshold_spin.value(), self.minpx_spin.value()
        )
        self._worker.finished.connect(self._on_recon_done)
        self._worker.error.connect(self._on_recon_error)
        self._worker.start()

    def _on_recon_done(self):
        self.recon_progress.setVisible(False)
        self.run_btn.setEnabled(True)
        if self.merge_check.isChecked():
            self.gb.merge(self.merge_angle_spin.value())
        n = self.gb.n_grains
        self.recon_info.setPlainText(
            f"Grains detected: {n}\n"
            f"Threshold: {self.threshold_spin.value()}°\n"
            f"Min pixels: {self.minpx_spin.value()}"
        )
        self._draw_grain_map()
        if self.ebsd_tab is not None:
            try:
                self.ebsd_tab.update_grains_summary(self.gb)
            except Exception:
                pass

    def _on_recon_error(self, msg: str):
        self.recon_progress.setVisible(False)
        self.run_btn.setEnabled(True)
        QMessageBox.critical(self, "Reconstruction error", msg)

    def _draw_grain_map(self):
        if not self.gb.is_ready:
            return
        grains = self.gb._grains
        data   = grains.data
        gids   = grains.grain_ids
        opts   = self._recon_map_opts
        ps     = opts.pt_size

        # --- Grain ID map ---
        fig1 = go.Figure(go.Scattergl(
            x=data["X"].values, y=data["Y"].values,
            mode="markers",
            marker=dict(
                color=gids.astype(float),
                colorscale="Turbo",
                size=ps,
                showscale=opts.show_colorbar,
                colorbar=dict(title="Grain ID") if opts.show_colorbar else None,
            ),
            hovertemplate="X=%{x:.1f}<br>Y=%{y:.1f}<br>GrainID=%{marker.color:.0f}<extra></extra>",
        ))
        fig1.update_layout(
            title="Grain ID map",
            xaxis_title="X (µm)", yaxis_title="Y (µm)",
            yaxis=dict(scaleanchor="x"),
        )
        self.grain_map_plt.show_figure(fig1)

        # --- Grain size map (area per grain) ---
        df_summary = self.gb.summary_df()
        if df_summary is not None and "area" in df_summary.columns:
            area_by_id = dict(zip(df_summary.index, df_summary["area"]))
            areas = np.array([area_by_id.get(g, 0.0) for g in gids], dtype=float)
        else:
            areas = gids.astype(float)

        vmin = opts.vmin if opts.vmin is not None else float(np.nanmin(areas[areas > 0]) if (areas > 0).any() else 0)
        vmax = opts.vmax if opts.vmax is not None else float(np.nanpercentile(areas, 99))

        fig2 = go.Figure(go.Scattergl(
            x=data["X"].values, y=data["Y"].values,
            mode="markers",
            marker=dict(
                color=areas,
                colorscale=opts.colorscale,
                cmin=vmin, cmax=vmax,
                size=ps,
                showscale=opts.show_colorbar,
                colorbar=dict(title="Area (px²)") if opts.show_colorbar else None,
            ),
            hovertemplate="X=%{x:.1f}<br>Y=%{y:.1f}<br>Area=%{marker.color:.1f}<extra></extra>",
        ))
        fig2.update_layout(
            title="Grain size map",
            xaxis_title="X (µm)", yaxis_title="Y (µm)",
            yaxis=dict(scaleanchor="x"),
        )
        self.size_map_plt.show_figure(fig2)

    # ==================================================================
    # Shape Parameters
    # ==================================================================

    def _build_shape_panel(self):
        panel = QWidget()
        layout = QHBoxLayout(panel)
        splitter = QSplitter(Qt.Horizontal)
        layout.addWidget(splitter)

        left = QWidget(); lv = QVBoxLayout(left)
        self.shape_load_btn = QPushButton("Compute shape parameters")
        self.shape_load_btn.clicked.connect(self._compute_shape)
        lv.addWidget(self.shape_load_btn)
        self.shape_table = QTableWidget()
        self.shape_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.shape_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        lv.addWidget(self.shape_table)
        splitter.addWidget(left)

        right = QWidget(); rv = QVBoxLayout(right)
        self.shape_vis = QTabWidget()
        self.hist_area_plt    = PlotlyWidget()
        self.hist_ar_plt      = PlotlyWidget()
        self.hist_gos_plt     = PlotlyWidget()
        self.scatter_shape_plt = PlotlyWidget()
        self.shape_vis.addTab(self.hist_area_plt,     "Area distribution")
        self.shape_vis.addTab(self.hist_ar_plt,       "Aspect ratio")
        self.shape_vis.addTab(self.hist_gos_plt,      "GOS distribution")
        self.shape_vis.addTab(self.scatter_shape_plt, "Size vs. shape")
        rv.addWidget(self.shape_vis)
        splitter.addWidget(right)
        splitter.setSizes([500, 600])
        return panel

    def _compute_shape(self):
        if not self.gb.is_ready:
            QMessageBox.warning(self, "No grains", "Run reconstruction first.")
            return
        df = self.gb.summary_df()
        if df is None:
            return

        # Table
        self.shape_table.clear()
        self.shape_table.setColumnCount(len(df.columns))
        self.shape_table.setHorizontalHeaderLabels(list(df.columns))
        self.shape_table.setRowCount(min(len(df), 2000))
        for row in range(min(len(df), 2000)):
            for col, cname in enumerate(df.columns):
                val = df.iloc[row][cname]
                txt = f"{val:.4f}" if isinstance(val, float) else str(val)
                self.shape_table.setItem(row, col, QTableWidgetItem(txt))

        # Area histogram
        fig_a = go.Figure(go.Histogram(
            x=df["area"], nbinsx=50,
            marker_color="#4e79a7",
            name="Area"
        ))
        fig_a.update_layout(title="Grain area distribution",
                            xaxis_title="Area (map units²)", yaxis_title="Count")
        self.hist_area_plt.show_figure(fig_a)

        # Aspect ratio histogram
        fig_ar = go.Figure(go.Histogram(
            x=df["aspect_ratio"], nbinsx=50,
            marker_color="#f28e2b"
        ))
        fig_ar.update_layout(title="Aspect ratio distribution",
                             xaxis_title="Aspect ratio (a/b)", yaxis_title="Count")
        self.hist_ar_plt.show_figure(fig_ar)

        # GOS histogram
        fig_gos = go.Figure(go.Histogram(
            x=df["gos_deg"], nbinsx=50,
            marker_color="#59a14f"
        ))
        fig_gos.update_layout(title="Grain Orientation Spread (GOS)",
                              xaxis_title="GOS (°)", yaxis_title="Count")
        self.hist_gos_plt.show_figure(fig_gos)

        # Scatter size vs shape
        fig_sc = go.Figure(go.Scatter(
            x=df["area"], y=df["aspect_ratio"],
            mode="markers",
            marker=dict(
                color=df["gos_deg"],
                colorscale="Plasma",
                size=5, opacity=0.7,
                colorbar=dict(title="GOS (°)"),
                showscale=True,
            ),
            hovertemplate="Area=%{x:.1f}<br>AR=%{y:.2f}<br>GOS=%{marker.color:.2f}°<extra></extra>",
        ))
        fig_sc.update_layout(title="Size vs. shape (colour = GOS)",
                             xaxis_title="Area", yaxis_title="Aspect ratio")
        self.scatter_shape_plt.show_figure(fig_sc)

    # ==================================================================
    # Orientation
    # ==================================================================

    def _build_orient_panel(self):
        panel = QWidget()
        layout = QHBoxLayout(panel)
        splitter = QSplitter(Qt.Horizontal)
        layout.addWidget(splitter)

        left = QWidget(); lv = QVBoxLayout(left)
        lv.addWidget(QLabel("Kernel Average Misorientation (KAM)"))
        kam_row = QHBoxLayout()
        kam_row.addWidget(QLabel("Max angle:"))
        self.kam_angle_spin = QDoubleSpinBox()
        self.kam_angle_spin.setRange(0.1, 20.0)
        self.kam_angle_spin.setSuffix("°")
        self.kam_angle_spin.setValue(5.0)
        kam_row.addWidget(self.kam_angle_spin)
        lv.addLayout(kam_row)

        self._kam_opts  = _PlotOptions("KAM options",  default_cmap="Hot", default_pt_size=2)
        self._grod_opts = _PlotOptions("GROD options", default_cmap="Inferno", default_pt_size=2)

        self.kam_btn = QPushButton("Compute KAM")
        self.kam_btn.clicked.connect(self._compute_kam)
        lv.addWidget(self.kam_btn)
        lv.addWidget(self._kam_opts)

        self.grod_btn = QPushButton("Compute GROD")
        self.grod_btn.clicked.connect(self._compute_grod)
        lv.addWidget(self.grod_btn)
        lv.addWidget(self._grod_opts)

        self.texture_label = QLabel("Texture index (J): —")
        self.texture_label.setWordWrap(True)
        lv.addWidget(self.texture_label)
        self.texture_btn = QPushButton("Compute texture index")
        self.texture_btn.clicked.connect(self._compute_texture)
        lv.addWidget(self.texture_btn)

        self.orient_progress = QProgressBar()
        self.orient_progress.setRange(0, 0)
        self.orient_progress.setVisible(False)
        lv.addWidget(self.orient_progress)
        lv.addStretch()
        splitter.addWidget(left)

        right = QWidget(); rv = QVBoxLayout(right)
        self.orient_vis = QTabWidget()
        self.kam_plt  = PlotlyWidget()
        self.grod_plt = PlotlyWidget()
        self.orient_vis.addTab(self.kam_plt,  "KAM map")
        self.orient_vis.addTab(self.grod_plt, "GROD map")
        rv.addWidget(self.orient_vis)
        splitter.addWidget(right)
        splitter.setSizes([320, 780])
        return panel

    def _compute_kam(self):
        if not self.gb.is_ready:
            QMessageBox.warning(self, "No grains", "Run reconstruction first.")
            return
        self.orient_progress.setVisible(True)
        max_a = self.kam_angle_spin.value()
        self._analysis_worker = _AnalysisWorker(lambda: self.gb.calc_kam(max_a))
        self._analysis_worker.finished.connect(self._on_kam_done)
        self._analysis_worker.error.connect(
            lambda m: (self.orient_progress.setVisible(False),
                       QMessageBox.critical(self, "KAM error", m)))
        self._analysis_worker.start()

    def _on_kam_done(self, kam):
        self.orient_progress.setVisible(False)
        data = self.gb._grains.data
        opts = self._kam_opts
        vmin = opts.vmin or 0.0
        vmax = opts.vmax or float(np.nanpercentile(kam, 99))
        cb = dict(title="KAM (°)") if opts.show_colorbar else None
        fig = go.Figure(go.Scattergl(
            x=data["X"].values, y=data["Y"].values,
            mode="markers",
            marker=dict(color=kam, colorscale=opts.colorscale,
                        cmin=vmin, cmax=vmax, size=opts.pt_size,
                        colorbar=cb, showscale=opts.show_colorbar),
        ))
        fig.update_layout(title="KAM map", xaxis_title="X", yaxis_title="Y",
                          yaxis=dict(scaleanchor="x"))
        self.kam_plt.show_figure(fig)
        self.orient_vis.setCurrentWidget(self.kam_plt)

    def _compute_grod(self):
        if not self.gb.is_ready:
            return
        self.orient_progress.setVisible(True)
        self._analysis_worker = _AnalysisWorker(self.gb.calc_grod)
        self._analysis_worker.finished.connect(self._on_grod_done)
        self._analysis_worker.error.connect(
            lambda m: (self.orient_progress.setVisible(False),
                       QMessageBox.critical(self, "GROD error", m)))
        self._analysis_worker.start()

    def _on_grod_done(self, grod):
        self.orient_progress.setVisible(False)
        data = self.gb._grains.data
        opts = self._grod_opts
        vmin = opts.vmin
        vmax = opts.vmax or float(np.nanpercentile(grod, 99))
        cb = dict(title="GROD (°)") if opts.show_colorbar else None
        fig = go.Figure(go.Scattergl(
            x=data["X"].values, y=data["Y"].values,
            mode="markers",
            marker=dict(color=grod, colorscale=opts.colorscale,
                        cmin=vmin, cmax=vmax, size=opts.pt_size,
                        colorbar=cb, showscale=opts.show_colorbar),
        ))
        fig.update_layout(title="GROD map", xaxis_title="X", yaxis_title="Y",
                          yaxis=dict(scaleanchor="x"))
        self.grod_plt.show_figure(fig)
        self.orient_vis.setCurrentWidget(self.grod_plt)

    def _compute_texture(self):
        if not self.gb.is_ready:
            return
        j = self.gb.texture_index()
        self.texture_label.setText(f"Texture index (J): {j:.3f}")

    # ==================================================================
    # Grain Boundaries
    # ==================================================================

    def _build_boundary_panel(self):
        panel = QWidget()
        layout = QHBoxLayout(panel)
        splitter = QSplitter(Qt.Horizontal)
        layout.addWidget(splitter)

        left = QWidget(); lv = QVBoxLayout(left)
        lagb_box = QGroupBox("Sub-grain boundaries")
        lf = QFormLayout(lagb_box)
        self.lagb_spin = QDoubleSpinBox()
        self.lagb_spin.setRange(0.1, 45.0)
        self.lagb_spin.setValue(15.0)
        self.lagb_spin.setSuffix("°")
        lf.addRow("LAGB threshold:", self.lagb_spin)
        lv.addWidget(lagb_box)

        csl_box = QGroupBox("CSL / Twinning (cubic)")
        clf = QFormLayout(csl_box)
        self.csl_spin = QDoubleSpinBox()
        self.csl_spin.setRange(0.5, 15.0)
        self.csl_spin.setValue(5.0)
        self.csl_spin.setSuffix("°")
        clf.addRow("CSL tolerance:", self.csl_spin)
        lv.addWidget(csl_box)

        self._bnd_opts = _PlotOptions("Boundary options", default_cmap="Hot",
                                       default_pt_size=2)
        lv.addWidget(self._bnd_opts)

        self.boundary_btn = QPushButton("Analyse boundaries")
        self.boundary_btn.clicked.connect(self._analyse_boundaries)
        lv.addWidget(self.boundary_btn)
        self.boundary_progress = QProgressBar()
        self.boundary_progress.setRange(0, 0)
        self.boundary_progress.setVisible(False)
        lv.addWidget(self.boundary_progress)
        self.boundary_info = QTextEdit()
        self.boundary_info.setReadOnly(True)
        lv.addWidget(self.boundary_info)
        lv.addStretch()
        splitter.addWidget(left)

        right = QWidget(); rv = QVBoxLayout(right)
        self.boundary_vis = QTabWidget()
        self.boundary_map_plt  = PlotlyWidget()
        self.mori_dist_plt     = PlotlyWidget()
        self.csl_map_plt       = PlotlyWidget()
        self.twist_tilt_plt    = PlotlyWidget()
        self.boundary_vis.addTab(self.boundary_map_plt, "Boundary map")
        self.boundary_vis.addTab(self.mori_dist_plt,    "Misorientation dist.")
        self.boundary_vis.addTab(self.csl_map_plt,      "CSL / twins")
        self.boundary_vis.addTab(self.twist_tilt_plt,   "Twist / tilt")
        rv.addWidget(self.boundary_vis)
        splitter.addWidget(right)
        splitter.setSizes([310, 790])
        return panel

    def _analyse_boundaries(self):
        if not self.gb.is_ready:
            QMessageBox.warning(self, "No grains", "Run reconstruction first.")
            return
        self.boundary_progress.setVisible(True)
        self.boundary_btn.setEnabled(False)
        lagb = self.lagb_spin.value()
        csl_thresh = self.csl_spin.value()

        def _run():
            gb = self.gb._grains.boundary
            mori_angles = gb.misorientation_angle
            lengths     = gb.segment_length
            n_seg       = gb.n_segments
            lagb_count  = int(np.sum(mori_angles < lagb))
            hagb_count  = n_seg - lagb_count
            total_len   = float(lengths.sum())
            lagb_len    = float(lengths[mori_angles < lagb].sum())
            sigma_vals  = gb.sigma_value(csl_thresh)
            twin_mask   = gb.is_twinning(3, csl_thresh)
            twin_len    = float(lengths[twin_mask].sum())
            tt          = gb.twist_tilt()
            centers, density = gb.angle_distribution(bins=36)
            return {
                "n_seg": n_seg, "lagb_count": lagb_count, "hagb_count": hagb_count,
                "total_len": total_len, "lagb_len": lagb_len, "twin_len": twin_len,
                "sigma_vals": sigma_vals, "twin_mask": twin_mask,
                "centers": centers, "density": density,
                "twist_angle": tt["twist_angle"], "tilt_angle": tt["tilt_angle"],
                "mori_angles": mori_angles, "lengths": lengths,
            }

        self._bnd_worker = _AnalysisWorker(_run)
        self._bnd_worker.finished.connect(self._on_boundary_done)
        self._bnd_worker.error.connect(
            lambda m: (self.boundary_progress.setVisible(False),
                       self.boundary_btn.setEnabled(True),
                       QMessageBox.critical(self, "Boundary error", m)))
        self._bnd_worker.start()

    def _on_boundary_done(self, result: dict):
        self.boundary_progress.setVisible(False)
        self.boundary_btn.setEnabled(True)
        opts = self._bnd_opts

        lines = [
            f"Total segments: {result['n_seg']}",
            f"LAGB (<{self.lagb_spin.value()}°): {result['lagb_count']}",
            f"HAGB: {result['hagb_count']}",
            f"Total length: {result['total_len']:.1f}",
            f"LAGB length: {result['lagb_len']:.1f} "
            f"({100*result['lagb_len']/max(result['total_len'],1e-9):.1f}%)",
            f"Sigma3 twin length: {result['twin_len']:.1f} "
            f"({100*result['twin_len']/max(result['total_len'],1e-9):.1f}%)",
        ]
        sigma_vals = result["sigma_vals"]
        for s in np.unique(sigma_vals[sigma_vals > 0]):
            n = int(np.sum(sigma_vals == s))
            lines.append(f"  Sigma{s}: {n} segments")
        self.boundary_info.setPlainText("\n".join(lines))

        # --- Boundary map: scatter at midpoints coloured by mori angle ---
        gb  = self.gb._grains.boundary
        mp  = gb.midpoint
        mori = result["mori_angles"]
        data = self.gb._grains.data
        vmin = opts.vmin or 0.0
        vmax = opts.vmax or float(np.nanpercentile(mori, 99))
        cb = dict(title="Misorientation (°)") if opts.show_colorbar else None
        fig_bnd = go.Figure()
        # Background: all EBSD points
        fig_bnd.add_trace(go.Scattergl(
            x=data["X"].values, y=data["Y"].values,
            mode="markers",
            marker=dict(color="#e0e0e0", size=1, opacity=0.3),
            name="EBSD points", showlegend=False,
        ))
        fig_bnd.add_trace(go.Scattergl(
            x=mp[:, 0], y=mp[:, 1],
            mode="markers",
            marker=dict(color=mori, colorscale=opts.colorscale,
                        cmin=vmin, cmax=vmax, size=opts.pt_size,
                        colorbar=cb, showscale=opts.show_colorbar),
            name="Boundaries",
            hovertemplate="Mori=%{marker.color:.1f}°<extra></extra>",
        ))
        fig_bnd.update_layout(title="Grain boundary map",
                              yaxis=dict(scaleanchor="x"))
        self.boundary_map_plt.show_figure(fig_bnd)

        # --- Misorientation distribution ---
        fig_mori = go.Figure()
        fig_mori.add_trace(go.Bar(
            x=result["centers"], y=result["density"],
            marker_color="#4e79a7", name="Length-weighted freq.",
        ))
        fig_mori.add_vline(x=self.lagb_spin.value(),
                           line_dash="dash", line_color="red",
                           annotation_text=f"LAGB ({self.lagb_spin.value()}°)")
        fig_mori.update_layout(title="Misorientation angle distribution",
                               xaxis_title="Misorientation angle (°)",
                               yaxis_title="Length-weighted frequency")
        self.mori_dist_plt.show_figure(fig_mori)

        # --- CSL map: scatter at midpoints, twin = red, HAGB = blue, rest = grey ---
        twin_mask   = result["twin_mask"]
        non_twin_hagb = ~twin_mask & gb.is_hagb(self.lagb_spin.value())
        fig_csl = go.Figure()
        fig_csl.add_trace(go.Scattergl(
            x=mp[~twin_mask & ~non_twin_hagb, 0],
            y=mp[~twin_mask & ~non_twin_hagb, 1],
            mode="markers", marker=dict(color="#cccccc", size=1, opacity=0.4),
            name="LAGB", showlegend=True,
        ))
        if non_twin_hagb.any():
            fig_csl.add_trace(go.Scattergl(
                x=mp[non_twin_hagb, 0], y=mp[non_twin_hagb, 1],
                mode="markers", marker=dict(color="#1f77b4", size=opts.pt_size),
                name="HAGB",
            ))
        if twin_mask.any():
            fig_csl.add_trace(go.Scattergl(
                x=mp[twin_mask, 0], y=mp[twin_mask, 1],
                mode="markers", marker=dict(color="red", size=opts.pt_size),
                name="Sigma3 twin",
            ))
        fig_csl.update_layout(title="CSL / twin boundaries",
                              yaxis=dict(scaleanchor="x"))
        self.csl_map_plt.show_figure(fig_csl)

        # --- Twist / tilt ---
        fig_tt = go.Figure(go.Scatter(
            x=result["twist_angle"], y=result["tilt_angle"],
            mode="markers",
            marker=dict(color="#9467bd", size=3, opacity=0.4),
        ))
        fig_tt.update_layout(title="Twist vs. Tilt decomposition",
                             xaxis_title="Twist angle (°)",
                             yaxis_title="Tilt angle (°)")
        self.twist_tilt_plt.show_figure(fig_tt)

    # ==================================================================
    # Triple Points
    # ==================================================================

    def _build_tp_panel(self):
        panel = QWidget()
        layout = QHBoxLayout(panel)
        splitter = QSplitter(Qt.Horizontal)
        layout.addWidget(splitter)

        left = QWidget(); lv = QVBoxLayout(left)
        self.tp_btn = QPushButton("Detect triple points")
        self.tp_btn.clicked.connect(self._detect_triple_points)
        lv.addWidget(self.tp_btn)
        self.tp_table = QTableWidget()
        self.tp_table.setEditTriggers(QAbstractItemView.NoEditTriggers)
        self.tp_table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeToContents)
        lv.addWidget(self.tp_table)
        lv.addStretch()
        splitter.addWidget(left)

        right = QWidget(); rv = QVBoxLayout(right)
        self.tp_plt = PlotlyWidget()
        rv.addWidget(self.tp_plt)
        splitter.addWidget(right)
        splitter.setSizes([450, 650])
        return panel

    def _detect_triple_points(self):
        if not self.gb.is_ready:
            QMessageBox.warning(self, "No grains", "Run reconstruction first.")
            return
        df = self.gb.triple_points_df()
        if df is None or len(df) == 0:
            QMessageBox.information(self, "Triple points", "No triple points detected.")
            return
        # Table
        self.tp_table.setColumnCount(len(df.columns))
        self.tp_table.setHorizontalHeaderLabels(list(df.columns))
        self.tp_table.setRowCount(min(len(df), 5000))
        for row in range(min(len(df), 5000)):
            for col, cn in enumerate(df.columns):
                v = df.iloc[row][cn]
                self.tp_table.setItem(row, col,
                    QTableWidgetItem(f"{v:.3f}" if isinstance(v, float) else str(v)))
        # Map
        grains = self.gb._grains
        data   = grains.data
        gids   = grains.grain_ids
        tp     = grains.triple_points
        fig = go.Figure()
        fig.add_trace(go.Scattergl(
            x=data["X"].values, y=data["Y"].values,
            mode="markers",
            marker=dict(color=gids.astype(float), colorscale="Turbo", size=1, opacity=0.4,
                        showscale=False),
            name="Grain IDs",
        ))
        if tp.n > 0:
            coords = tp.coordinates
            fig.add_trace(go.Scatter(
                x=coords[:, 0], y=coords[:, 1],
                mode="markers",
                marker=dict(color="red", size=6, symbol="x"),
                name="Triple points",
            ))
        fig.update_layout(title=f"Triple points ({tp.n})",
                          yaxis=dict(scaleanchor="x"))
        self.tp_plt.show_figure(fig)

    # ==================================================================
    # Export
    # ==================================================================

    def _build_export_panel(self):
        panel = QWidget()
        layout = QVBoxLayout(panel)
        layout.addWidget(QLabel("Export grain analysis results:"))
        for label, slot in (
            ("Export grain summary (CSV)",      self._export_csv),
            ("Export Neper orientations (.ori)", self._export_neper_ori),
            ("Export Neper tessellation (.tess)",self._export_neper_tess),
            ("Export EBSD + grain IDs (CSV)",   self._export_ctf),
        ):
            row = QWidget()
            rl = QHBoxLayout(row)
            btn = QPushButton(label)
            btn.clicked.connect(slot)
            rl.addWidget(btn)
            rl.addStretch()
            layout.addWidget(row)
        layout.addStretch()
        return panel

    def _guard(self) -> bool:
        if not self.gb.is_ready:
            QMessageBox.warning(self, "No grains", "Run reconstruction first.")
            return False
        return True

    def _export_csv(self):
        if not self._guard(): return
        path, _ = QFileDialog.getSaveFileName(self, "Save CSV", "", "CSV (*.csv)")
        if path:
            try:
                self.gb.export_csv(path)
                QMessageBox.information(self, "Exported", f"Saved to {path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e))

    def _export_neper_ori(self):
        if not self._guard(): return
        path, _ = QFileDialog.getSaveFileName(self, "Save .ori", "", "Neper ori (*.ori)")
        if path:
            try:
                self.gb.export_neper_ori(path)
                QMessageBox.information(self, "Exported", f"Saved to {path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e))

    def _export_neper_tess(self):
        if not self._guard(): return
        path, _ = QFileDialog.getSaveFileName(self, "Save .tess", "", "Neper tess (*.tess)")
        if path:
            try:
                self.gb.export_neper_tess(path)
                QMessageBox.information(self, "Exported", f"Saved to {path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e))

    def _export_ctf(self):
        if not self._guard(): return
        path, _ = QFileDialog.getSaveFileName(self, "Save EBSD+GrainID CSV", "", "CSV (*.csv)")
        if path:
            try:
                self.gb.export_ctf(path)
                QMessageBox.information(self, "Exported", f"Saved to {path}")
            except Exception as e:
                QMessageBox.critical(self, "Error", str(e))
