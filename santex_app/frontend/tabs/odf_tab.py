"""
ODFTab — ODF and Pole Figure analysis tab (Plotly visualisation).

Sub-tabs (visualisation)
------------------------
  ODF Sections    : phi2 / phi1 / sigma Euler-section contour plots
  Pole Figure     : stereographic pole figure with KDE contours
  IPF Map         : inverse pole figure coloured scatter map
  Components      : horizontal bar chart of texture component fractions
  Properties      : texture index, entropy, ODF max
  Export          : POPLA / MTEX-txt / VPSC
"""

from __future__ import annotations
import numpy as np

from PyQt5.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QLabel,
    QDoubleSpinBox, QSpinBox, QPushButton, QSplitter,
    QComboBox, QFormLayout, QTabWidget, QProgressBar,
    QCheckBox, QFileDialog, QMessageBox, QTextEdit,
    QLineEdit, QScrollArea,
)
from PyQt5.QtCore import Qt, QThread, pyqtSignal
from PyQt5.QtGui import QFont

from frontend.widgets.plotly_widget import PlotlyWidget, COLORSCALES
from frontend.tabs.ebsd_tab import _PlotOptions

import plotly.graph_objects as go
from plotly.subplots import make_subplots


# ---------------------------------------------------------------------------
# Worker thread
# ---------------------------------------------------------------------------

class _ODFWorker(QThread):
    finished = pyqtSignal(str)
    error    = pyqtSignal(str)

    def __init__(self, task: str, fn):
        super().__init__()
        self.task = task
        self.fn   = fn

    def run(self):
        try:
            self.fn()
            self.finished.emit(self.task)
        except Exception as e:
            self.error.emit(f"{self.task}: {e}")


# ---------------------------------------------------------------------------
# Symmetry display list — populated lazily
# ---------------------------------------------------------------------------

_SYM_DISPLAY_LIST: list[tuple[str, str]] = []

_SECTION_TYPES    = ["phi2 sections", "phi1 sections", "sigma sections"]
_SECTION_DEFAULTS = {
    "phi2 sections":  [0, 15, 30, 45, 60, 75, 90],
    "phi1 sections":  [0, 15, 30, 45, 60, 75, 90],
    "sigma sections": [0, 15, 30, 45, 60, 75, 90],
}


# ---------------------------------------------------------------------------
# ODFTab
# ---------------------------------------------------------------------------

class ODFTab(QWidget):
    def __init__(self, odf_backend, ebsd_backend, parent=None):
        super().__init__(parent)
        self.ob = odf_backend
        self.eb = ebsd_backend
        self._worker = None
        self._build_ui()

    # ------------------------------------------------------------------

    def _build_ui(self):
        # Build symmetry list once
        global _SYM_DISPLAY_LIST
        try:
            _SYM_DISPLAY_LIST = self.ob.symmetry_display_list()
        except Exception:
            _SYM_DISPLAY_LIST = [
                ("Cubic (m-3m)  [Oh]  —  24 ops",        "m-3m"),
                ("Hexagonal (6/mmm)  [D6h]  —  12 ops",  "6/mmm"),
                ("Orthorhombic (mmm)  [D2h]  —  4 ops",  "mmm"),
                ("Tetragonal (4/mmm)  [D4h]  —  8 ops",  "4/mmm"),
                ("Trigonal (-3m)  [D3d]  —  6 ops",      "-3m"),
                ("Monoclinic (2/m)  [C2h]  —  2 ops",    "2/m"),
                ("Triclinic (-1)  [Ci]  —  1 ops",       "-1"),
            ]

        root = QHBoxLayout(self)
        splitter = QSplitter(Qt.Horizontal)
        root.addWidget(splitter)

        # ---- Left: scrollable controls ----
        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFixedWidth(370)
        left = QWidget()
        lv = QVBoxLayout(left)
        lv.setContentsMargins(4, 4, 4, 4)
        scroll.setWidget(left)

        # ODF computation
        comp_box = QGroupBox("ODF Computation")
        comp_form = QFormLayout(comp_box)

        self.phase_spin = QSpinBox()
        self.phase_spin.setRange(0, 20)
        self.phase_spin.setValue(1)
        comp_form.addRow("Phase index:", self.phase_spin)

        self.sym_combo = QComboBox()
        self.sym_combo.setToolTip(
            "All 11 Laue groups supported.\n"
            "Format: Name  [Schoenflies]  — n proper rotations."
        )
        for label, _key in _SYM_DISPLAY_LIST:
            self.sym_combo.addItem(label)
        _default_keys = [k for _, k in _SYM_DISPLAY_LIST]
        _mmm_idx = _default_keys.index("mmm") if "mmm" in _default_keys else 0
        self.sym_combo.setCurrentIndex(_mmm_idx)
        comp_form.addRow("Crystal symmetry:", self.sym_combo)

        self.ssym_combo = QComboBox()
        self.ssym_combo.addItems(["none", "orthorhombic"])
        comp_form.addRow("Sample symmetry:", self.ssym_combo)

        self.hw_spin = QDoubleSpinBox()
        self.hw_spin.setRange(1.0, 45.0)
        self.hw_spin.setValue(10.0)
        self.hw_spin.setSuffix("°")
        self.hw_spin.setSingleStep(1.0)
        comp_form.addRow("Kernel halfwidth:", self.hw_spin)

        self.res_spin = QDoubleSpinBox()
        self.res_spin.setRange(2.0, 15.0)
        self.res_spin.setValue(5.0)
        self.res_spin.setSuffix("°")
        comp_form.addRow("Grid resolution:", self.res_spin)

        self.max_ori_spin = QSpinBox()
        self.max_ori_spin.setRange(500, 50000)
        self.max_ori_spin.setValue(5000)
        self.max_ori_spin.setSingleStep(500)
        comp_form.addRow("Max orientations:", self.max_ori_spin)

        self.compute_btn = QPushButton("Compute ODF")
        self.compute_btn.clicked.connect(self._start_compute_odf)
        comp_form.addRow(self.compute_btn)

        self.progress = QProgressBar()
        self.progress.setRange(0, 0)
        self.progress.setVisible(False)
        comp_form.addRow(self.progress)
        lv.addWidget(comp_box)

        # ODF section controls
        sec_box = QGroupBox("ODF Sections")
        sec_form = QFormLayout(sec_box)
        self.sec_type_combo = QComboBox()
        self.sec_type_combo.addItems(_SECTION_TYPES)
        sec_form.addRow("Section type:", self.sec_type_combo)
        self.sec_vals_edit = QLineEdit("0, 15, 30, 45, 60, 75, 90")
        sec_form.addRow("Section values (°):", self.sec_vals_edit)
        self._sec_opts = _PlotOptions("Section plot options", default_cmap="RdYlBu_r",
                                       default_pt_size=2)
        sec_form.addRow(self._sec_opts)
        self.plot_sections_btn = QPushButton("Plot Sections")
        self.plot_sections_btn.clicked.connect(self._plot_sections)
        sec_form.addRow(self.plot_sections_btn)
        lv.addWidget(sec_box)

        # Pole figure
        pf_box = QGroupBox("Pole Figure")
        pf_form = QFormLayout(pf_box)
        hkl_row = QHBoxLayout()
        self.pf_h = QSpinBox(); self.pf_h.setRange(-10, 10); self.pf_h.setValue(1)
        self.pf_k = QSpinBox(); self.pf_k.setRange(-10, 10); self.pf_k.setValue(0)
        self.pf_l = QSpinBox(); self.pf_l.setRange(-10, 10); self.pf_l.setValue(0)
        for w in (self.pf_h, self.pf_k, self.pf_l):
            w.setMaximumWidth(55); hkl_row.addWidget(w)
        pf_form.addRow("h k l:", hkl_row)
        self.pf_bw_spin = QDoubleSpinBox()
        self.pf_bw_spin.setRange(1.0, 30.0)
        self.pf_bw_spin.setValue(7.5)
        self.pf_bw_spin.setSuffix("°")
        pf_form.addRow("KDE bandwidth:", self.pf_bw_spin)
        self._pf_opts = _PlotOptions("PF plot options", default_cmap="RdYlBu_r",
                                      default_pt_size=2)
        pf_form.addRow(self._pf_opts)
        self.pf_btn = QPushButton("Plot Pole Figure")
        self.pf_btn.clicked.connect(self._plot_pf)
        pf_form.addRow(self.pf_btn)
        lv.addWidget(pf_box)

        # IPF
        ipf_box = QGroupBox("Inverse Pole Figure")
        ipf_form = QFormLayout(ipf_box)
        self.ipf_dir_combo = QComboBox()
        self.ipf_dir_combo.addItems(["Z (ND)", "X (RD)", "Y (TD)"])
        ipf_form.addRow("Sample direction:", self.ipf_dir_combo)
        self.ipf_pt_size = QSpinBox()
        self.ipf_pt_size.setRange(1, 10)
        self.ipf_pt_size.setValue(2)
        ipf_form.addRow("Point size:", self.ipf_pt_size)
        self.ipf_btn = QPushButton("Plot IPF Map")
        self.ipf_btn.clicked.connect(self._plot_ipf)
        ipf_form.addRow(self.ipf_btn)
        lv.addWidget(ipf_box)

        # Component / export
        comp_ana_btn = QPushButton("Analyse Components")
        comp_ana_btn.clicked.connect(self._analyse_components)
        lv.addWidget(comp_ana_btn)

        prop_btn = QPushButton("Compute Properties (J, S)")
        prop_btn.clicked.connect(self._compute_properties)
        lv.addWidget(prop_btn)

        export_box = QGroupBox("Export ODF")
        ev = QVBoxLayout(export_box)
        for label, slot in (
            ("Export POPLA (.epf)",      self._export_popla),
            ("Export MTEX txt",          self._export_mtex),
            ("Export VPSC orientations", self._export_vpsc),
        ):
            b = QPushButton(label)
            b.clicked.connect(slot)
            ev.addWidget(b)
        lv.addWidget(export_box)
        lv.addStretch()
        splitter.addWidget(scroll)

        # ---- Right: visualisation ----
        right = QWidget()
        rv = QVBoxLayout(right)
        self.vis_tabs = QTabWidget()

        self.sections_plt = PlotlyWidget()
        self.vis_tabs.addTab(self.sections_plt, "ODF Sections")

        self.pf_plt = PlotlyWidget()
        self.vis_tabs.addTab(self.pf_plt, "Pole Figure")

        self.ipf_plt = PlotlyWidget()
        self.vis_tabs.addTab(self.ipf_plt, "IPF Map")

        self.comp_plt = PlotlyWidget()
        self.vis_tabs.addTab(self.comp_plt, "Components")

        self.prop_text = QTextEdit()
        self.prop_text.setReadOnly(True)
        mono = QFont("Courier New", 9)
        mono.setStyleHint(QFont.Monospace)
        self.prop_text.setFont(mono)
        self.vis_tabs.addTab(self.prop_text, "Properties")

        rv.addWidget(self.vis_tabs)
        splitter.addWidget(right)
        splitter.setSizes([370, 730])

    # ------------------------------------------------------------------
    # ODF Computation
    # ------------------------------------------------------------------

    def _start_compute_odf(self):
        if not self.eb.is_loaded:
            QMessageBox.warning(self, "No EBSD data", "Load an EBSD file first.")
            return
        idx     = self.sym_combo.currentIndex()
        sym_str = _SYM_DISPLAY_LIST[idx][1] if idx < len(_SYM_DISPLAY_LIST) else "mmm"
        phase   = self.phase_spin.value()
        hw      = self.hw_spin.value()
        max_ori = self.max_ori_spin.value()
        data    = self.eb.data

        def _fn():
            self.ob.compute_odf(data, phase=phase,
                                crystal_symmetry=sym_str,
                                halfwidth_deg=hw,
                                max_orientations=max_ori)

        self._run_worker("ODF", _fn)

    def _run_worker(self, task: str, fn):
        self.compute_btn.setEnabled(False)
        self.progress.setVisible(True)
        self._worker = _ODFWorker(task, fn)
        self._worker.finished.connect(self._on_worker_done)
        self._worker.error.connect(self._on_worker_error)
        self._worker.start()

    def _on_worker_done(self, task: str):
        self.progress.setVisible(False)
        self.compute_btn.setEnabled(True)
        if task == "ODF":
            QMessageBox.information(self, "ODF ready",
                                    "ODF computation complete.\n"
                                    "Use the buttons to plot sections, PF, or IPF.")

    def _on_worker_error(self, msg: str):
        self.progress.setVisible(False)
        self.compute_btn.setEnabled(True)
        QMessageBox.critical(self, "Error", msg)

    # ------------------------------------------------------------------
    # ODF Sections
    # ------------------------------------------------------------------

    def _parse_section_values(self) -> list[float]:
        try:
            return [float(v.strip()) for v in self.sec_vals_edit.text().split(",") if v.strip()]
        except ValueError:
            return [0, 15, 30, 45, 60, 75, 90]

    def _plot_sections(self):
        if not self.ob.is_ready:
            QMessageBox.warning(self, "No ODF", "Compute ODF first.")
            return

        sec_type = self.sec_type_combo.currentText()
        vals     = self._parse_section_values()
        res      = self.res_spin.value()
        opts     = self._sec_opts

        if sec_type == "phi2 sections":
            sections = self.ob.get_phi2_sections(vals, res)
            xlabel, ylabel, fixed_label = "phi1 (deg)", "PHI (deg)", "phi2"
        elif sec_type == "phi1 sections":
            sections = self.ob.get_phi1_sections(vals, res)
            xlabel, ylabel, fixed_label = "phi2 (deg)", "PHI (deg)", "phi1"
        else:
            sections = self.ob.get_sigma_sections(vals, res)
            xlabel, ylabel, fixed_label = "phi2 (deg)", "PHI (deg)", "sigma"

        if not sections:
            return

        n     = len(sections)
        ncols = min(4, n)
        nrows = (n + ncols - 1) // ncols

        all_odf = np.concatenate([v[2].ravel() for v in sections.values()])
        vmax_data = float(np.percentile(all_odf, 99))
        vmax = opts.vmax if opts.vmax else vmax_data
        vmin = opts.vmin if opts.vmin is not None else 0.0

        # Subplot titles
        subtitles = [f"{fixed_label} = {sv:.0f}°" for sv in sections.keys()]

        fig = make_subplots(rows=nrows, cols=ncols, subplot_titles=subtitles,
                            shared_xaxes=False, shared_yaxes=False,
                            horizontal_spacing=0.06, vertical_spacing=0.10)

        # Build shared colorbar (only on last trace)
        all_keys = list(sections.keys())
        for i, (phi_val, (p1_2d, P_2d, odf_2d)) in enumerate(sections.items()):
            row = i // ncols + 1
            col = i % ncols  + 1
            is_last = (i == len(sections) - 1)

            fig.add_trace(go.Contour(
                z=odf_2d.T,
                x=p1_2d[0, :],
                y=P_2d[:, 0],
                colorscale=opts.colorscale,
                zmin=vmin, zmax=vmax,
                contours=dict(start=vmin, end=vmax,
                              size=(vmax - vmin) / 12 if vmax > vmin else 1),
                colorbar=dict(title="ODF (m.r.d.)", x=1.02) if (is_last and opts.show_colorbar) else None,
                showscale=is_last and opts.show_colorbar,
                line_smoothing=0.8,
                line_width=0.5,
                name=f"{fixed_label}={phi_val:.0f}°",
                showlegend=False,
            ), row=row, col=col)

            # Axis labels only on outer subplots
            fig.update_xaxes(title_text=xlabel if row == nrows else "", row=row, col=col)
            fig.update_yaxes(title_text=ylabel if col == 1 else "", row=row, col=col)

        fig.update_layout(
            title=f"ODF {sec_type}  (halfwidth={self.ob.halfwidth_deg:.0f}°, "
                  f"symmetry={self.ob.crystal_symmetry})",
            height=max(400, nrows * 220),
        )
        self.sections_plt.show_figure(fig)
        self.vis_tabs.setCurrentWidget(self.sections_plt)

    # ------------------------------------------------------------------
    # Pole Figure
    # ------------------------------------------------------------------

    def _plot_pf(self):
        if not self.ob.is_ready:
            QMessageBox.warning(self, "No ODF", "Compute ODF first.")
            return

        hkl    = (self.pf_h.value(), self.pf_k.value(), self.pf_l.value())
        bw     = self.pf_bw_spin.value()
        ssym   = self.ssym_combo.currentText()
        opts   = self._pf_opts

        result = self.ob.get_pole_figure_kde(hkl, ssym, bw)
        if result is None:
            return
        xy_grid, density, xy_scatter = result

        vmin = opts.vmin if opts.vmin is not None else 0.0
        vmax = opts.vmax if opts.vmax else float(np.nanpercentile(density, 99))
        cb   = dict(title="m.r.d.") if opts.show_colorbar else None

        fig = go.Figure()
        # Reference circle
        th = np.linspace(0, 2 * np.pi, 300)
        fig.add_shape(type="circle", xref="x", yref="y",
                      x0=-1, y0=-1, x1=1, y1=1,
                      line=dict(color="black", width=1.5))
        # Cross-hairs
        for x0, x1, y0, y1 in [(-1, 1, 0, 0), (0, 0, -1, 1)]:
            fig.add_shape(type="line", x0=x0, x1=x1, y0=y0, y1=y1,
                          line=dict(color="gray", width=0.5, dash="dot"))

        # KDE contour
        try:
            from scipy.interpolate import griddata
            xi = np.linspace(-1, 1, 150)
            yi = np.linspace(-1, 1, 150)
            xi_g, yi_g = np.meshgrid(xi, yi)
            zi = griddata(xy_grid, density, (xi_g, yi_g), method="linear",
                          fill_value=0.0)
            mask = xi_g ** 2 + yi_g ** 2 > 1
            zi[mask] = np.nan

            fig.add_trace(go.Contour(
                z=zi, x=xi, y=yi,
                colorscale=opts.colorscale,
                zmin=vmin, zmax=vmax,
                contours=dict(start=vmin, end=vmax,
                              size=(vmax - vmin) / 10 if vmax > vmin else 0.1),
                colorbar=cb,
                showscale=opts.show_colorbar,
                line_smoothing=0.8,
            ))
        except Exception:
            # Fallback: scatter density
            fig.add_trace(go.Scattergl(
                x=xy_grid[:, 0], y=xy_grid[:, 1],
                mode="markers",
                marker=dict(color=density, colorscale=opts.colorscale,
                            size=4, colorbar=cb, showscale=opts.show_colorbar),
            ))

        # Raw pole scatter
        fig.add_trace(go.Scattergl(
            x=xy_scatter[:, 0], y=xy_scatter[:, 1],
            mode="markers",
            marker=dict(color="rgba(80,80,80,0.15)", size=opts.pt_size),
            name="raw PF", showlegend=False,
        ))

        # N / S / E / W labels
        for text, x, y, ax_, ay_ in [
            ("N", 0, 1.1, 0, 15), ("S", 0, -1.1, 0, -15),
            ("E", 1.1, 0, 15, 0), ("W", -1.1, 0, -15, 0),
        ]:
            fig.add_annotation(x=x, y=y, text=text, showarrow=False,
                               font=dict(size=11))

        fig.update_layout(
            title=f"Pole Figure ({hkl[0]} {hkl[1]} {hkl[2]})  —  {self.ob.crystal_symmetry}",
            xaxis=dict(range=[-1.15, 1.15], showgrid=False, zeroline=False,
                       showticklabels=False, scaleanchor="y"),
            yaxis=dict(range=[-1.15, 1.15], showgrid=False, zeroline=False,
                       showticklabels=False),
        )
        self.pf_plt.show_figure(fig)
        self.vis_tabs.setCurrentWidget(self.pf_plt)

    # ------------------------------------------------------------------
    # IPF Map
    # ------------------------------------------------------------------

    def _plot_ipf(self):
        if not self.ob.is_ready:
            QMessageBox.warning(self, "No ODF", "Compute ODF first.")
            return
        if self.eb.data is None:
            return

        dir_map    = {"Z (ND)": (0, 0, 1), "X (RD)": (1, 0, 0), "Y (TD)": (0, 1, 0)}
        sample_dir = dir_map[self.ipf_dir_combo.currentText()]
        rgb        = self.ob.get_ipf_colors(sample_dir)
        if rgb is None:
            return

        data  = self.eb.data
        phase = self.ob._phase_id
        mask  = data["Phase"] == phase if phase is not None else data["Phase"] > 0
        x_vals = data.loc[mask, "X"].values
        y_vals = data.loc[mask, "Y"].values
        ps     = self.ipf_pt_size.value()

        # Convert (N,3) float RGB → hex strings for Plotly
        rgb_u8  = np.clip(rgb * 255, 0, 255).astype(np.uint8)
        colors  = [f"rgb({r},{g},{b})" for r, g, b in rgb_u8]

        fig = make_subplots(rows=1, cols=2, column_widths=[0.78, 0.22],
                            subplot_titles=["IPF map", "IPF key"])

        fig.add_trace(go.Scattergl(
            x=x_vals, y=y_vals,
            mode="markers",
            marker=dict(color=colors, size=ps),
            hovertemplate="X=%{x:.1f}<br>Y=%{y:.1f}<extra></extra>",
        ), row=1, col=1)
        fig.update_xaxes(title_text="X (µm)", row=1, col=1)
        fig.update_yaxes(title_text="Y (µm)", scaleanchor="x", row=1, col=1)

        # IPF colour triangle
        tri = self.ob.get_ipf_triangle_colors(n_pts=3000)
        if tri is not None:
            xy_tri, rgb_tri = tri
            rgb_u8_tri  = np.clip(rgb_tri * 255, 0, 255).astype(np.uint8)
            colors_tri  = [f"rgb({r},{g},{b})" for r, g, b in rgb_u8_tri]
            fig.add_trace(go.Scatter(
                x=xy_tri[:, 0], y=xy_tri[:, 1],
                mode="markers",
                marker=dict(color=colors_tri, size=4),
                showlegend=False,
            ), row=1, col=2)
            # Corner labels
            for text, tx, ty in [("001", 0, 0), ("011", 0.5, 0.5), ("111", 0.35, 0.35)]:
                fig.add_annotation(x=tx, y=ty, text=text, showarrow=False,
                                   font=dict(size=9), row=1, col=2)
        fig.update_xaxes(showticklabels=False, showgrid=False, row=1, col=2)
        fig.update_yaxes(showticklabels=False, showgrid=False, row=1, col=2)

        fig.update_layout(
            title=f"IPF map — {self.ipf_dir_combo.currentText()}  ({self.ob.crystal_symmetry})",
            showlegend=False,
        )
        self.ipf_plt.show_figure(fig)
        self.vis_tabs.setCurrentWidget(self.ipf_plt)

    # ------------------------------------------------------------------
    # Texture Components
    # ------------------------------------------------------------------

    def _analyse_components(self):
        if not self.ob.is_ready:
            QMessageBox.warning(self, "No ODF", "Compute ODF first.")
            return

        rows  = self.ob.component_analysis(max_angle_deg=15.0)
        if not rows:
            return

        names = [r["name"]        for r in rows]
        fracs = [r["fraction"]    for r in rows]
        descs = [r["description"] for r in rows]
        pcts  = [f * 100          for f in fracs]

        # Horizontal bar with description annotations
        colors = [f"hsl({int(i/len(names)*300)},70%,55%)" for i in range(len(names))]
        fig = go.Figure(go.Bar(
            x=pcts, y=names,
            orientation="h",
            marker_color=colors,
            text=[f"  {d}" for d in descs],
            textposition="outside",
            hovertemplate="<b>%{y}</b><br>%{x:.1f}%<extra></extra>",
        ))
        fig.update_layout(
            title=f"Texture components (within 15° — {self.ob.crystal_symmetry})",
            xaxis_title="Volume fraction (%)",
            xaxis=dict(range=[0, max(pcts) * 1.6 + 1 if pcts else 20]),
            margin=dict(l=120, r=200),
        )
        self.comp_plt.show_figure(fig)
        self.vis_tabs.setCurrentWidget(self.comp_plt)

    # ------------------------------------------------------------------
    # ODF Properties
    # ------------------------------------------------------------------

    def _compute_properties(self):
        if not self.ob.is_ready:
            QMessageBox.warning(self, "No ODF", "Compute ODF first.")
            return

        self.progress.setVisible(True)
        self.compute_btn.setEnabled(False)

        def _do():
            res = 10.0
            J    = self.ob.texture_index(res)
            S    = self.ob.entropy(res)
            fmax = self.ob.odf_max(5.0)
            lines = [
                "ODF Properties",
                "=" * 40,
                f"Crystal symmetry : {self.ob.crystal_symmetry}",
                f"Kernel halfwidth  : {self.ob.halfwidth_deg:.1f}°",
                f"Phase             : {self.ob._phase_id}",
                f"No. orientations  : {self.ob._odf._N}",
                "",
                f"Texture index (J) : {J:.4f}  m.r.d.²",
                "  (J=1 random, J>1 textured)",
                "",
                f"Entropy (S)       : {S:.4f}",
                "  (S=0 random, S<0 textured)",
                "",
                f"ODF max           : {fmax:.2f}  m.r.d.",
                "",
                "Volume fractions (within 15°):",
            ]
            for r in self.ob.component_analysis(15.0):
                lines.append(f"  {r['name']:<16s}  {r['fraction']*100:5.1f}%  {r['description']}")
            self.prop_text.setPlainText("\n".join(lines))

        self._worker = _ODFWorker("Properties", _do)
        self._worker.finished.connect(lambda _: (
            self.progress.setVisible(False),
            self.compute_btn.setEnabled(True),
            self.vis_tabs.setCurrentWidget(self.prop_text),
        ))
        self._worker.error.connect(self._on_worker_error)
        self._worker.start()

    # ------------------------------------------------------------------
    # Export
    # ------------------------------------------------------------------

    def _export_popla(self):
        self._export("POPLA files (*.epf);;All files (*)", self.ob.export_popla)

    def _export_mtex(self):
        self._export("Text files (*.txt);;All files (*)", self.ob.export_mtex_txt)

    def _export_vpsc(self):
        self._export("Text files (*.txt);;All files (*)", self.ob.export_vpsc)

    def _export(self, filt: str, fn):
        if not self.ob.is_ready:
            QMessageBox.warning(self, "No ODF", "Compute ODF first.")
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export ODF", "", filt)
        if not path:
            return
        try:
            fn(path)
            QMessageBox.information(self, "Export", f"Saved to {path}")
        except Exception as e:
            QMessageBox.critical(self, "Export error", str(e))
