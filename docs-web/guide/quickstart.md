# Quick Start

## Workflow 1 — Single mineral anisotropy

The fastest way to see SAnTex in action:

1. Launch the app (`python santex_app/main.py`).
2. Open the **Material Database** tab.
3. Click **Forsterite** in the mineral table.
4. Leave P = 0 GPa, T = 300 K and click **Compute anisotropy**.
5. A *Vp* stereonet appears instantly on the right.

You will see the Voigt matrix, scalar anisotropy metrics (AVS, Vp max/min), and an interactive Plotly stereonet you can zoom, hover, and export as SVG.

## Workflow 2 — Textured rock from EBSD

```
Open EBSD file  →  Filter  →  Map phases  →  VRH average  →  Plot
```

1. **File → Open EBSD file…** and select a `.ctf` scan.
2. The **EBSD** tab opens automatically with a phase map.
3. *(Optional)* Set a rectangular ROI, apply MAD ≤ 0.8, set downsample factor.
4. In the **VRH Averaging** section, map each phase number to its mineral name.
5. Set pressure, temperature, and averaging method (Voigt / Reuss / Hill).
6. Click **Compute VRH** — the textured elastic tensor and *Vp* stereonet appear.
7. Switch to **Grains** → **Reconstruct grains** for per-grain shape and misorientation statistics.
8. Switch to **ODF & PF** → set crystal symmetry → **Compute ODF** → plot φ₂ sections and pole figures.

## Workflow 3 — Modal rock (no EBSD)

1. Open the **Modal Rock** tab.
2. Add rows: phase name + volume fraction. Fractions must sum to 100 %.
3. Set P (GPa) and T (K).
4. Click **Compute Modal Rock** — Voigt-averaged tensor + *Vp* stereonet.
5. Switch to the **Hashin–Shtrikman** sub-tab for isotropic upper/lower/VRH bounds.

## Scripting (no GUI)

```python
from santex.material import Material
from santex.anisotropy import Anisotropy

mat = Material()

# Forsterite at 2 GPa, 1000 °C
cij = mat.voigt_high_PT("Forsterite", PRESSURE=2, TEMP=1273)
rho = mat.load_density("Forsterite", pressure=2, temperature=1273)

anis = Anisotropy(cij * 1e9, rho)   # units: Pa, kg/m³
anis.anisotropy_values()             # print scalar metrics
anis.plot()                          # matplotlib figure
```

::: tip
All core computations (`santex/` library) work without Qt. Use them in Jupyter notebooks or scripts independently of the GUI.
:::
