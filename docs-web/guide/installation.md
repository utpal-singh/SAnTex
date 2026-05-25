# Installation

## Requirements

| Package | Version | Purpose |
|---------|---------|---------|
| Python | ≥ 3.10 | Runtime |
| numpy | ≥ 1.23 | Array math |
| scipy | ≥ 1.9 | Linear algebra, interpolation |
| pandas | ≥ 1.5 | EBSD data frames |
| orix | ≥ 0.12 | Crystal orientations |
| plotly | ≥ 5.11 | 2-D interactive plots |
| PyQt5 | ≥ 5.15 | GUI framework |
| PyQtWebEngine | ≥ 5.15 | Plotly chart renderer |
| pyvista | ≥ 0.40 | 3-D velocity surface |
| pyvistaqt | ≥ 0.11 | PyVista Qt integration |
| vtk | ≥ 9.1 | 3-D engine |
| scikit-learn | ≥ 1.1 | ML hole-filling |
| tabulate | ≥ 0.9 | Table formatting |

## Conda (recommended)

```bash
conda create -n santex python=3.10
conda activate santex
pip install -r requirements.txt
pip install PyQtWebEngine
```

::: warning Windows — PyQt5 / PyQtWebEngine DLL conflict
conda-forge's PyQt5 ships `Qt5Core_conda.dll` (with a `_conda` suffix) while pip's
PyQtWebEngine expects `Qt5Core.dll`.  Fix this by replacing the entire PyQt5 stack
with pip versions:

```bash
pip install --force-reinstall PyQt5 PyQtWebEngine
```
:::

## pip (any virtual environment)

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install numpy scipy pandas orix plotly \
            PyQt5 PyQtWebEngine \
            pyvista pyvistaqt vtk \
            scikit-learn tabulate tqdm joblib
```

## Run the application

```bash
cd santex_app
python main.py
```

Or from the repo root using conda:

```bash
conda run -n santex --no-capture-output python santex_app/main.py
```

## Verify the installation

```python
# Quick smoke test — no GUI
import santex
from santex.material import Material
m = Material()
print(m.available_phases())       # prints the full mineral table

from santex.anisotropy import Anisotropy
import numpy as np
cij = m.get_voigt_matrix("Forsterite") * 1e9   # GPa → Pa
rho = m.load_density("Forsterite")
anis = Anisotropy(cij, rho)
vals = anis.anisotropy_values()
print(f"Vp anisotropy: {vals['p_wave_anisotropy_percent']:.2f} %")
```

Expected output:
```
Vp anisotropy: 24.68 %
```
