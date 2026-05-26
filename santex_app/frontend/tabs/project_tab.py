"""
ProjectTab — multi-dataset EBSD project browser.

UI layout
---------
  Left   : dataset list (QTreeWidget with checkboxes, cache badge, phase info)
  Right  : detail panel for selected dataset + batch operations
  Bottom : memory settings (max RAM slider) + project file controls

Signals
-------
  open_dataset(dataset_id: str)   — emitted when user requests a dataset
                                    to open in the EBSD analysis tab
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import QColor, QFont, QIcon
from PyQt5.QtWidgets import (
    QAbstractItemView,
    QAction,
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMenu,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSlider,
    QSpinBox,
    QSplitter,
    QTextEdit,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
    QProgressDialog,
    QApplication,
)


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

_CACHE_STYLE = "color:#27ae60; font-weight:bold;"   # green "●" when in cache
_DISK_STYLE  = "color:#95a5a6;"                      # grey  "○" when on disk


def _humanise_pixels(n: int) -> str:
    if n >= 1_000_000:
        return f"{n/1e6:.2f} M px"
    if n >= 1_000:
        return f"{n/1e3:.0f} K px"
    return f"{n} px"


def _phase_sym_labels(meta) -> str:
    """One-liner like 'Forsterite (mmm), Enstatite (mmm)'."""
    from santex.ebsd.ipf_coloring import spacegroup_to_laue
    from backend.ebsd_backend import _LAUE_KEY_LABEL

    if meta.phases_df is None or meta.phases_df.empty:
        return "—"
    parts = []
    for idx in meta.phases_df.index:
        row   = meta.phases_df.loc[idx]
        pname = str(row.get("phase", f"Phase{idx}")).strip()
        try:
            laue_key = spacegroup_to_laue(row.get("symmetry", 11))
            hm       = _LAUE_KEY_LABEL.get(laue_key, laue_key)
        except Exception:
            hm = "?"
        parts.append(f"{pname} ({hm})")
    return ";  ".join(parts)


# ─────────────────────────────────────────────────────────────────────────────
# DatasetTreeWidget  — the main list
# ─────────────────────────────────────────────────────────────────────────────

class DatasetTreeWidget(QTreeWidget):
    """QTreeWidget that shows EBSD datasets with checkboxes + cache badge.

    Columns
    -------
    0  ☑ Name (editable display name)
    1  Pixels
    2  Step (µm)
    3  Phases (symmetry)
    4  Cache
    """

    # Emitted on double-click
    open_requested = pyqtSignal(str)   # dataset_id

    def __init__(self, project_tab=None, parent=None):
        super().__init__(parent)
        self._project_tab = project_tab  # back-reference set by ProjectTab after construction
        self.setColumnCount(5)
        self.setHeaderLabels(["Dataset", "Pixels", "Step", "Phase(s)", "RAM"])
        self.setAlternatingRowColors(True)
        self.setSelectionMode(QAbstractItemView.SingleSelection)
        self.setUniformRowHeights(True)
        self.setContextMenuPolicy(Qt.CustomContextMenu)
        self.customContextMenuRequested.connect(self._ctx_menu)
        self.itemDoubleClicked.connect(self._on_double_click)

        hdr = self.header()
        hdr.setSectionResizeMode(0, QHeaderView.Stretch)
        hdr.setSectionResizeMode(1, QHeaderView.ResizeToContents)
        hdr.setSectionResizeMode(2, QHeaderView.ResizeToContents)
        hdr.setSectionResizeMode(3, QHeaderView.ResizeToContents)
        hdr.setSectionResizeMode(4, QHeaderView.Fixed)
        self.setColumnWidth(4, 44)

        # id → QTreeWidgetItem
        self._items: dict[str, QTreeWidgetItem] = {}

    # ── public API ─────────────────────────────────────────────────────

    def add_dataset(self, meta) -> None:
        """Add or update a row for *meta*."""
        item = self._items.get(meta.id)
        if item is None:
            item = QTreeWidgetItem(self)
            item.setData(0, Qt.UserRole, meta.id)
            item.setFlags(
                Qt.ItemIsEnabled | Qt.ItemIsSelectable | Qt.ItemIsUserCheckable
            )
            self._items[meta.id] = item
        self._refresh_item(item, meta)

    def remove_dataset(self, dataset_id: str) -> None:
        item = self._items.pop(dataset_id, None)
        if item:
            idx = self.indexOfTopLevelItem(item)
            if idx >= 0:
                self.takeTopLevelItem(idx)

    def refresh_dataset(self, meta) -> None:
        item = self._items.get(meta.id)
        if item:
            self._refresh_item(item, meta)

    def refresh_all(self, datasets: list) -> None:
        for meta in datasets:
            item = self._items.get(meta.id)
            if item:
                self._refresh_item(item, meta)

    def selected_id(self) -> Optional[str]:
        items = self.selectedItems()
        if items:
            return items[0].data(0, Qt.UserRole)
        return None

    def all_ids(self) -> list[str]:
        return list(self._items.keys())

    # ── internals ──────────────────────────────────────────────────────

    def _refresh_item(self, item: QTreeWidgetItem, meta) -> None:
        item.setCheckState(0, Qt.Checked if meta.active else Qt.Unchecked)
        item.setText(0, meta.name)
        item.setText(1, _humanise_pixels(meta.n_pixels))
        item.setText(2, f"{meta.step_um:.3g} µm" if meta.step_um else "—")
        item.setText(3, _phase_sym_labels(meta))
        item.setText(4, "●" if meta.in_cache else "○")
        item.setForeground(4, QColor("#27ae60") if meta.in_cache else QColor("#95a5a6"))
        item.setToolTip(0, str(meta.path))
        item.setTextAlignment(1, Qt.AlignRight | Qt.AlignVCenter)
        item.setTextAlignment(2, Qt.AlignRight | Qt.AlignVCenter)
        item.setTextAlignment(4, Qt.AlignCenter)

    def _on_double_click(self, item: QTreeWidgetItem, _col: int) -> None:
        ds_id = item.data(0, Qt.UserRole)
        if ds_id:
            self.open_requested.emit(ds_id)

    def _ctx_menu(self, pos) -> None:
        item = self.itemAt(pos)
        if item is None:
            return
        ds_id = item.data(0, Qt.UserRole)
        menu  = QMenu(self)

        open_act   = menu.addAction("Open in EBSD Tab")
        rename_act = menu.addAction("Rename…")
        menu.addSeparator()
        evict_act  = menu.addAction("Release from RAM")
        menu.addSeparator()
        remove_act = menu.addAction("Remove from project")

        act = menu.exec_(self.viewport().mapToGlobal(pos))
        if act == open_act:
            self.open_requested.emit(ds_id)
        elif act == rename_act:
            self._rename(item, ds_id)
        elif act == evict_act:
            if self._project_tab:
                self._project_tab._evict(ds_id)
        elif act == remove_act:
            if self._project_tab:
                self._project_tab._remove_dataset(ds_id)

    def _rename(self, item: QTreeWidgetItem, ds_id: str) -> None:
        dlg = QDialog(self)
        dlg.setWindowTitle("Rename dataset")
        vb  = QVBoxLayout(dlg)
        le  = QLineEdit(item.text(0))
        vb.addWidget(le)
        bb = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        bb.accepted.connect(dlg.accept)
        bb.rejected.connect(dlg.reject)
        vb.addWidget(bb)
        if dlg.exec_() == QDialog.Accepted and le.text().strip():
            item.setText(0, le.text().strip())
            if self._project_tab:
                self._project_tab._on_renamed(ds_id, le.text().strip())


# ─────────────────────────────────────────────────────────────────────────────
# ProjectTab
# ─────────────────────────────────────────────────────────────────────────────

class ProjectTab(QWidget):
    """Multi-dataset EBSD project browser.

    Signals
    -------
    open_dataset(dataset_id: str)
        Emitted when the user wants to open a dataset in the EBSD tab.
    project_changed()
        Emitted whenever the dataset list or active flags change.
    """

    open_dataset    = pyqtSignal(str)
    project_changed = pyqtSignal()

    def __init__(self, project, ebsd_backend, parent=None):
        super().__init__(parent)
        self._proj  = project        # EBSDProject
        self._eb    = ebsd_backend   # EBSDBackend
        self._build_ui()

    # ── UI construction ────────────────────────────────────────────────

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(6)

        # ── Top bar: project controls ──────────────────────────────────
        top_bar = QHBoxLayout()
        top_bar.setSpacing(6)

        self._project_lbl = QLabel("Project: <b>Untitled</b>")
        top_bar.addWidget(self._project_lbl, 1)

        self._new_btn  = QPushButton("New")
        self._open_btn = QPushButton("Open…")
        self._save_btn = QPushButton("Save")
        self._saveas_btn = QPushButton("Save As…")
        for btn, slot in [
            (self._new_btn,    self._new_project),
            (self._open_btn,   self._open_project),
            (self._save_btn,   self._save_project),
            (self._saveas_btn, self._saveas_project),
        ]:
            btn.setFixedWidth(76)
            btn.clicked.connect(slot)
            top_bar.addWidget(btn)

        root.addLayout(top_bar)

        # ── Splitter: dataset list | detail panel ──────────────────────
        splitter = QSplitter(Qt.Horizontal)

        # LEFT — dataset list
        left = QWidget()
        lv   = QVBoxLayout(left)
        lv.setContentsMargins(0, 0, 0, 0)
        lv.setSpacing(4)

        self._tree = DatasetTreeWidget(project_tab=self, parent=self)
        self._tree.open_requested.connect(self._on_open_requested)
        self._tree.itemSelectionChanged.connect(self._on_selection_changed)
        self._tree.itemChanged.connect(self._on_item_changed)
        lv.addWidget(self._tree, 1)

        # Add / Remove buttons
        btn_row = QHBoxLayout()
        self._add_btn      = QPushButton("+ Add files…")
        self._add_dir_btn  = QPushButton("+ Scan folder…")
        self._remove_btn   = QPushButton("Remove selected")
        self._remove_btn.setEnabled(False)
        for btn, slot in [
            (self._add_btn,     self._add_files),
            (self._add_dir_btn, self._add_folder),
            (self._remove_btn,  lambda: self._remove_dataset(self._tree.selected_id())),
        ]:
            btn.clicked.connect(slot)
            btn_row.addWidget(btn)
        lv.addLayout(btn_row)

        # Memory indicator
        self._mem_lbl = QLabel()
        self._mem_lbl.setStyleSheet("color: #555; font-size: 11px;")
        lv.addWidget(self._mem_lbl)

        splitter.addWidget(left)

        # RIGHT — detail + batch panel (scroll area)
        right_scroll = QScrollArea()
        right_scroll.setWidgetResizable(True)
        right_scroll.setMinimumWidth(300)
        right_w = QWidget()
        right_v = QVBoxLayout(right_w)
        right_v.setContentsMargins(4, 4, 4, 4)
        right_v.setSpacing(8)
        right_scroll.setWidget(right_w)

        # ── Detail panel ───────────────────────────────────────────────
        detail_box = QGroupBox("Selected dataset")
        detail_form = QFormLayout(detail_box)
        detail_form.setSpacing(4)

        self._det_name   = QLabel("—")
        self._det_path   = QLabel("—")
        self._det_path.setWordWrap(True)
        self._det_path.setStyleSheet("font-size: 10px; color:#555;")
        self._det_grid   = QLabel("—")
        self._det_step   = QLabel("—")
        self._det_npix   = QLabel("—")
        self._det_phases = QLabel("—")
        self._det_phases.setWordWrap(True)
        self._det_cache  = QLabel("—")

        for lbl, widget in [
            ("Name:",    self._det_name),
            ("Path:",    self._det_path),
            ("Grid:",    self._det_grid),
            ("Step:",    self._det_step),
            ("Pixels:",  self._det_npix),
            ("Phase(s):", self._det_phases),
            ("In RAM:",  self._det_cache),
        ]:
            detail_form.addRow(lbl, widget)

        self._open_in_ebsd_btn = QPushButton("Open in EBSD Tab  →")
        self._open_in_ebsd_btn.setEnabled(False)
        self._open_in_ebsd_btn.clicked.connect(
            lambda: self._on_open_requested(self._tree.selected_id())
        )
        detail_form.addRow(self._open_in_ebsd_btn)
        right_v.addWidget(detail_box)

        # ── Memory settings ────────────────────────────────────────────
        mem_box  = QGroupBox("Memory settings")
        mem_form = QFormLayout(mem_box)
        mem_form.setSpacing(4)

        self._ram_spin = QSpinBox()
        self._ram_spin.setRange(1, 32)
        self._ram_spin.setValue(self._proj.max_in_ram)
        self._ram_spin.setSuffix(" dataset(s) in RAM")
        self._ram_spin.setToolTip(
            "Maximum number of full EBSD DataFrames held in memory simultaneously.\n"
            "Increasing this makes dataset switching faster at the cost of RAM.\n"
            "The least-recently-used entry is evicted when capacity is exceeded."
        )
        self._ram_spin.valueChanged.connect(self._on_ram_changed)
        mem_form.addRow("Max in RAM:", self._ram_spin)

        self._evict_btn = QPushButton("Release all from RAM")
        self._evict_btn.clicked.connect(self._evict_all)
        mem_form.addRow(self._evict_btn)
        right_v.addWidget(mem_box)

        # ── Batch operations ───────────────────────────────────────────
        batch_box  = QGroupBox("Batch operations  (apply to all ☑ datasets)")
        batch_form = QFormLayout(batch_box)
        batch_form.setSpacing(4)
        batch_form.addRow(QLabel(
            "<small><i>These load each active dataset in turn, apply the\n"
            "operation, and save the result back.</i></small>"
        ))

        self._batch_mad_spin = QDoubleSpinBox()
        self._batch_mad_spin.setRange(0.0, 10.0)
        self._batch_mad_spin.setValue(0.7)
        self._batch_mad_spin.setDecimals(2)
        self._batch_mad_spin.setSuffix("°")
        batch_form.addRow("MAD threshold:", self._batch_mad_spin)

        self._batch_mad_btn = QPushButton("Apply MAD filter to all ☑")
        self._batch_mad_btn.clicked.connect(self._batch_mad)
        batch_form.addRow(self._batch_mad_btn)

        batch_form.addRow(QLabel("<hr>"))

        self._batch_export_btn = QPushButton("Export all ☑ as CSV…")
        self._batch_export_btn.clicked.connect(self._batch_export_csv)
        batch_form.addRow(self._batch_export_btn)

        right_v.addWidget(batch_box)
        right_v.addStretch()

        splitter.addWidget(right_scroll)
        splitter.setSizes([600, 340])
        root.addWidget(splitter, 1)

        # Initial state
        self._refresh_all()
        self._update_mem_label()

    # ── project-level controls ──────────────────────────────────────────

    def _new_project(self):
        from backend.ebsd_project import EBSDProject
        if self._proj.all_datasets():
            reply = QMessageBox.question(
                self, "New project",
                "Discard the current project and start fresh?",
                QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
            )
            if reply != QMessageBox.Yes:
                return
        new_proj = EBSDProject(max_in_ram=self._proj.max_in_ram)
        self._replace_project(new_proj)

    def _open_project(self):
        from backend.ebsd_project import EBSDProject
        path, _ = QFileDialog.getOpenFileName(
            self, "Open project", "",
            f"SAnTex project (*{EBSDProject.FILE_EXTENSION});;All files (*)",
        )
        if not path:
            return
        try:
            proj = EBSDProject.load_project(path)
            self._replace_project(proj)
        except Exception as e:
            QMessageBox.critical(self, "Open project", f"Failed to open project:\n{e}")

    def _save_project(self):
        if not self._proj.is_saved:
            self._saveas_project()
            return
        try:
            self._proj.save(self._proj.project_path)
            self._update_title()
        except Exception as e:
            QMessageBox.critical(self, "Save project", f"Failed to save:\n{e}")

    def _saveas_project(self):
        from backend.ebsd_project import EBSDProject
        path, _ = QFileDialog.getSaveFileName(
            self, "Save project as", "",
            f"SAnTex project (*{EBSDProject.FILE_EXTENSION});;All files (*)",
        )
        if not path:
            return
        if not path.endswith(EBSDProject.FILE_EXTENSION):
            path += EBSDProject.FILE_EXTENSION
        try:
            self._proj.save(path)
            self._update_title()
        except Exception as e:
            QMessageBox.critical(self, "Save project", f"Failed to save:\n{e}")

    def _replace_project(self, new_proj) -> None:
        """Swap the active project object everywhere."""
        self._proj = new_proj
        self._eb._project = new_proj   # keep backend in sync (informal reference)
        # Rebuild list widget
        self._tree.clear()
        self._tree._items.clear()
        for meta in new_proj.all_datasets():
            self._tree.add_dataset(meta)
        self._ram_spin.setValue(new_proj.max_in_ram)
        self._update_title()
        self._update_mem_label()
        self._clear_detail()
        self.project_changed.emit()

    def _update_title(self):
        pp = self._proj.project_path
        name = Path(pp).stem if pp else "Untitled"
        self._project_lbl.setText(f"Project: <b>{name}</b>")

    # ── adding / removing datasets ──────────────────────────────────────

    def _add_files(self):
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Add EBSD files", "",
            "EBSD files (*.ctf *.ang *.osc);;All files (*)",
        )
        if not paths:
            return
        for path in paths:
            try:
                meta = self._proj.add_file(path)
                self._tree.add_dataset(meta)
            except Exception as e:
                QMessageBox.warning(self, "Add file", f"Could not add {path}:\n{e}")
        self._update_mem_label()
        self.project_changed.emit()

    def _add_folder(self):
        folder = QFileDialog.getExistingDirectory(self, "Scan folder for EBSD files")
        if not folder:
            return
        found = list(Path(folder).rglob("*.ctf")) + list(Path(folder).rglob("*.ang"))
        if not found:
            QMessageBox.information(self, "Scan folder", "No CTF/ANG files found.")
            return
        added = 0
        for path in sorted(found):
            try:
                meta = self._proj.add_file(str(path))
                self._tree.add_dataset(meta)
                added += 1
            except Exception:
                pass
        QMessageBox.information(
            self, "Scan folder",
            f"Added {added} file(s) from\n{folder}",
        )
        self._update_mem_label()
        self.project_changed.emit()

    def _remove_dataset(self, dataset_id: Optional[str]):
        if not dataset_id:
            return
        meta = self._proj.get_meta(dataset_id)
        name = meta.name if meta else dataset_id
        reply = QMessageBox.question(
            self, "Remove dataset",
            f"Remove '{name}' from the project?\n(The file on disk is not deleted.)",
            QMessageBox.Yes | QMessageBox.No, QMessageBox.No,
        )
        if reply != QMessageBox.Yes:
            return
        self._proj.remove(dataset_id)
        self._tree.remove_dataset(dataset_id)
        self._clear_detail()
        self._remove_btn.setEnabled(False)
        self._update_mem_label()
        self.project_changed.emit()

    def _on_renamed(self, dataset_id: str, new_name: str) -> None:
        self._proj.rename(dataset_id, new_name)
        self.project_changed.emit()

    # ── cache / RAM controls ────────────────────────────────────────────

    def _evict(self, dataset_id: str) -> None:
        self._proj.evict(dataset_id)
        meta = self._proj.get_meta(dataset_id)
        if meta:
            self._tree.refresh_dataset(meta)
        self._update_mem_label()

    def _evict_all(self) -> None:
        self._proj.evict_all()
        self._tree.refresh_all(self._proj.all_datasets())
        self._update_mem_label()

    def _on_ram_changed(self, value: int) -> None:
        self._proj.max_in_ram = value
        self._tree.refresh_all(self._proj.all_datasets())
        self._update_mem_label()

    # ── selection / detail panel ────────────────────────────────────────

    def _on_selection_changed(self):
        ds_id = self._tree.selected_id()
        self._remove_btn.setEnabled(ds_id is not None)
        self._open_in_ebsd_btn.setEnabled(ds_id is not None)
        if ds_id:
            meta = self._proj.get_meta(ds_id)
            if meta:
                self._show_detail(meta)
        else:
            self._clear_detail()

    def _on_item_changed(self, item: QTreeWidgetItem, col: int):
        if col == 0:
            ds_id  = item.data(0, Qt.UserRole)
            active = item.checkState(0) == Qt.Checked
            self._proj.set_active(ds_id, active)
            self.project_changed.emit()

    def _show_detail(self, meta) -> None:
        self._det_name.setText(f"<b>{meta.name}</b>")
        self._det_path.setText(str(meta.path))
        self._det_grid.setText(f"{meta.xcells} × {meta.ycells}")
        self._det_step.setText(f"{meta.step_um:.4g} µm" if meta.step_um else "—")
        self._det_npix.setText(_humanise_pixels(meta.n_pixels))
        self._det_phases.setText(_phase_sym_labels(meta))
        self._det_cache.setText(
            '<span style="color:#27ae60; font-weight:bold;">Yes — in RAM</span>'
            if meta.in_cache else
            '<span style="color:#888;">No — on disk</span>'
        )

    def _clear_detail(self) -> None:
        for w in (self._det_name, self._det_path, self._det_grid,
                  self._det_step, self._det_npix, self._det_phases, self._det_cache):
            w.setText("—")

    # ── open in EBSD tab ───────────────────────────────────────────────

    def _on_open_requested(self, dataset_id: Optional[str]) -> None:
        if not dataset_id:
            return
        self.open_dataset.emit(dataset_id)

    # ── batch operations ────────────────────────────────────────────────

    def _batch_mad(self):
        active = self._proj.active_datasets()
        if not active:
            QMessageBox.information(self, "Batch", "No active (☑) datasets.")
            return
        threshold = self._batch_mad_spin.value()

        prog = QProgressDialog(
            f"Applying MAD ≤ {threshold}° …", "Cancel", 0, len(active), self
        )
        prog.setWindowModality(Qt.WindowModal)
        prog.show()

        from santex.ebsd import EBSD as _EBSD
        for i, meta in enumerate(active):
            if prog.wasCanceled():
                break
            prog.setLabelText(f"Processing {meta.name}…")
            prog.setValue(i)
            QApplication.processEvents()
            try:
                df = self._proj.get_data(meta.id)
                if df is None:
                    continue
                ebsd = _EBSD(meta.path)
                filtered = ebsd.filter_MAD(df, threshold)
                # Put filtered copy back into cache (replaces raw)
                self._proj._cache.put(meta.id, filtered)
            except Exception as e:
                QMessageBox.warning(self, "Batch MAD", f"{meta.name}: {e}")

        prog.setValue(len(active))
        prog.close()
        QMessageBox.information(self, "Batch", f"MAD filter applied to {len(active)} dataset(s).")

    def _batch_export_csv(self):
        active = self._proj.active_datasets()
        if not active:
            QMessageBox.information(self, "Batch export", "No active (☑) datasets.")
            return
        folder = QFileDialog.getExistingDirectory(self, "Export CSV files to folder")
        if not folder:
            return

        prog = QProgressDialog("Exporting…", "Cancel", 0, len(active), self)
        prog.setWindowModality(Qt.WindowModal)
        prog.show()

        for i, meta in enumerate(active):
            if prog.wasCanceled():
                break
            prog.setLabelText(f"Exporting {meta.name}…")
            prog.setValue(i)
            QApplication.processEvents()
            try:
                df = self._proj.get_data(meta.id)
                if df is None:
                    continue
                out = Path(folder) / f"{meta.name}.csv"
                df.to_csv(str(out), index=False)
            except Exception as e:
                QMessageBox.warning(self, "Export", f"{meta.name}: {e}")

        prog.setValue(len(active))
        prog.close()
        QMessageBox.information(self, "Batch export", f"Exported {len(active)} file(s) to\n{folder}")

    # ── helpers ────────────────────────────────────────────────────────

    def _refresh_all(self):
        for meta in self._proj.all_datasets():
            self._tree.add_dataset(meta)

    def _update_mem_label(self):
        n = self._proj.cached_count
        cap = self._proj.max_in_ram
        total = len(self._proj)
        self._mem_lbl.setText(
            f"  {total} dataset(s) in project  ·  "
            f"{n}/{cap} slot(s) used  (RAM budget)"
        )

    # ── called by MainWindow after a dataset is loaded via project ─────

    def refresh_cache_badges(self):
        """Refresh the ● / ○ cache badges for all rows."""
        self._tree.refresh_all(self._proj.all_datasets())
        self._update_mem_label()
        # Also refresh the detail panel if something is selected
        ds_id = self._tree.selected_id()
        if ds_id:
            meta = self._proj.get_meta(ds_id)
            if meta:
                self._show_detail(meta)
