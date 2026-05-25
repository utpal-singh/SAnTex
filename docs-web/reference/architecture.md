# Architecture

## Directory layout

```
SAnTex/
├── santex/                        ← Core library (GUI-independent)
│   ├── material/material.py
│   ├── ebsd/
│   │   ├── ebsd.py
│   │   ├── calcGrainBoundary.py
│   │   ├── odf.py
│   │   └── ...
│   ├── grains/
│   │   ├── grain2d.py
│   │   └── symmetry.py
│   ├── anisotropy/anisotropy.py
│   ├── odf/odf.py
│   ├── isotropy/isotropy.py
│   └── tensor/tensor.py
│
├── santex_app/                    ← Desktop GUI
│   ├── main.py                    ← Entry point
│   ├── frontend/
│   │   ├── main_window.py
│   │   ├── tabs/
│   │   │   ├── material_tab.py
│   │   │   ├── anisotropy_tab.py
│   │   │   ├── ebsd_tab.py        ← also defines _PlotOptions
│   │   │   ├── modal_rock_tab.py
│   │   │   ├── grains_tab.py
│   │   │   └── odf_tab.py
│   │   └── widgets/
│   │       ├── plotly_widget.py   ← QWebEngineView wrapper
│   │       └── pyvista_widget.py  ← BackgroundPlotter wrapper
│   └── backend/
│       ├── material_backend.py
│       ├── ebsd_backend.py
│       ├── grains_backend.py
│       ├── anisotropy_backend.py
│       └── odf_backend.py
│
└── docs-web/                      ← This website
```

## Layered architecture

```
┌──────────────────────────────────────────────────┐
│  GUI Widgets  (PyQt5, Plotly, PyVista)           │
│    ↕ signal/slot                                 │
│  Tab classes  (frontend/tabs/*.py)               │
│    ↕ method calls                                │
│  Backend classes  (backend/*.py)                 │
│    ↕ plain Python calls                          │
│  Core library  (santex/)  — zero Qt dependency  │
└──────────────────────────────────────────────────┘
```

**Core library** (`santex/`) has no Qt dependency and can be used from Jupyter notebooks or scripts.

**Backends** are thin adapters that translate between the core library's types and what the GUI tabs need. They do no plotting.

**Tabs** own the Qt widget layout and call backends. Long computations are offloaded to `QThread` workers so the GUI stays responsive.

**Widgets** are reusable rendering components shared across tabs.

## Worker thread pattern

Every long computation (grain reconstruction, ODF, VRH average, velocity surface) follows this pattern:

```python
class _Worker(QThread):
    finished = pyqtSignal(object)
    error    = pyqtSignal(str)

    def __init__(self, backend, **kwargs):
        super().__init__()
        self.backend = backend
        self.kwargs  = kwargs

    def run(self):
        try:
            result = self.backend.compute_something(**self.kwargs)
            self.finished.emit(result)
        except Exception as e:
            self.error.emit(str(e))

# In the tab:
self._worker = _Worker(self.backend, phase=phase, pressure=p)
self._worker.finished.connect(self._on_result)
self._worker.error.connect(self._on_error)
self._worker.start()
```

Keeping a reference (`self._worker`) prevents premature garbage collection.

## PlotlyWidget rendering pipeline

```
go.Figure  →  pio.to_html(include_plotlyjs=False)
           →  inject JS shim (fixes Qt 5.15 CSS insertRule bug)
           →  inject <script src="file:///TMPDIR/plotly.min.js">
           →  write to TMPDIR/plot_<id>.html
           →  QWebEngineView.load(QUrl.fromLocalFile(...))
```

`plotly.min.js` is **copied** from the Plotly package into the same temp directory as the HTML files at startup. This ensures both share the same `file://` origin, satisfying Qt WebEngine's same-origin security policy.

## Signal / slot connections map

| Source | Signal | Slot | Description |
|--------|--------|------|-------------|
| File menu | `triggered` | `MainWindow._open_file` | Open EBSD dialog |
| Tab selector | `currentChanged` | — | Switches visible tab |
| _Worker | `finished` | `Tab._on_done` | Delivers result to UI thread |
| _Worker | `error` | `Tab._on_error` | Shows error label |
| Compute button | `clicked` | `Tab._compute` | Starts worker |
