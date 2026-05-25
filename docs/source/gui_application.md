# SAnTex — Seismic Anisotropy Toolbox

**SAnTex** (Seismic Anisotropy from Texture) is an open-source Python application and library for computing the full elastic tensor of rocks from crystallographic texture data and calculating the resulting seismic wave velocities and anisotropy. It combines a material database of single-crystal elastic constants, EBSD (Electron Backscatter Diffraction) data processing, grain reconstruction, Orientation Distribution Function (ODF) analysis, and interactive 2-D / 3-D visualisation — all in a single desktop GUI.

> **Citation:** Singh *et al.* (2025), *Journal of Open Source Software*. Licensed under GNU GPL v3.

## Table of Contents

1. [Overview](#overview)
2. [Installation](#installation)
3. [Quick Start](#quick-start)
4. [Application Layout](#application-layout)
5. [Tabs — User Guide](#tabs-user-guide)
   - [Material Database](#material-database-tab)
   - [Anisotropy](#anisotropy-tab)
   - [EBSD](#ebsd-tab)
   - [Modal Rock](#modal-rock-tab)
   - [Grains](#grains-tab)
   - [ODF & PF](#odf-pf-tab)
6. [Python API Reference](#python-api-reference)
   - [Material Backend](#material-material-py)
   - [EBSD Backend](#ebsd-ebsd-py)
   - [Grains Backend](#grains-grains-py)
   - [Anisotropy Backend](#anisotropy-anisotropy-py)
   - [ODF Backend](#odf-odf-py)
   - [Isotropy](#isotropy-isotropy-py)
   - [Tensor Utilities](#tensor-tensor-py)
   - [Crystal Symmetry](#grains-symmetry-py)
7. [Architecture](#architecture)
8. [File Formats](#file-formats)
9. [Keyboard Shortcuts](#keyboard-shortcuts)
10. [Troubleshooting](#troubleshooting)

## Overview

### What SAnTex does

| Step | Input | Output |
|------|-------|--------|
| **1. Single-crystal database** | Phase name, P, T | Voigt stiffness matrix *C*ᵢⱼ (GPa) |
| **2. Texture averaging** | EBSD Euler angles + *C*ᵢⱼ | Polycrystal tensor (Voigt / Reuss / Hill) |
| **3. Velocity computation** | Polycrystal *C*ᵢⱼ, density ρ | *Vp*, *Vs*₁, *Vs*₂ surfaces (km/s) |
| **4. Anisotropy metrics** | Velocity surfaces | AVS, P-wave anisotropy, Vp/Vs ratios |
| **5. Modal rock** | Volume fractions + P, T | Rock-average tensor + isotropic HS bounds |

### Key features

- **Material database** — 40+ rock-forming minerals with pressure and temperature derivatives (Birch-Murnaghan EOS)
- **EBSD processing** — load `.ctf` / `.ang` / `.osc` files, filter by MAD, define spatial ROI, downsample
- **Grain reconstruction** — misorientation-based segmentation, shape analysis (ellipse, Feret, convex hull), KAM, GROD, CSL boundary detection
- **ODF / Pole figures** — kernel-density ODF on SO(3), all 11 Laue-group symmetries, KDE pole figures, IPF maps
- **Interactive plots** — Plotly (2-D stereonets, contour ODF sections, bar charts) + PyVista (3-D velocity surface)
- **Export** — POPLA `.epf`, MTEX-compatible text, VPSC orientation files, grain CSV, Neper TESS/ORI, `.ctf`

## Installation

### Requirements

| Package | Version |
|---------|---------|
| Python | ≥ 3.10 |
| numpy | ≥ 1.23 |
| scipy | ≥ 1.9 |
| pandas | ≥ 1.5 |
| orix | ≥ 0.12 |
| plotly | ≥ 5.11 |
| PyQt5 | ≥ 5.15 |
| PyQtWebEngine | ≥ 5.15 |
| pyvista | ≥ 0.40 |
| pyvistaqt | ≥ 0.11 |
| vtk | ≥ 9.1 |
| scikit-learn | ≥ 1.1 |
| tabulate | ≥ 0.9 |

### Conda (recommended)

```bash
conda create -n santex python=3.10
conda activate santex
pip install -r requirements.txt
pip install PyQtWebEngine
```

> **Windows note:** conda-forge's PyQt5 uses `_conda`-suffixed Qt DLLs that are incompatible with PyQtWebEngine from pip. Fix by replacing the entire PyQt5 stack with pip:
> ```bash
> pip install --force-reinstall PyQt5 PyQtWebEngine
> ```

### pip (any virtual environment)

```bash
pip install numpy scipy pandas orix plotly PyQt5 PyQtWebEngine \
            pyvista pyvistaqt vtk scikit-learn tabulate tqdm joblib
```

### Run the application

```bash
cd santex_app
python main.py
```

## Quick Start

### Computing anisotropy for a single mineral

1. Open the **Material Database** tab.
2. Click a mineral (e.g. *Forsterite*) in the table.
3. Set **Pressure** and **Temperature**.
4. Click **Compute anisotropy**.
5. The Voigt matrix, scalar metrics, and a *Vp* stereonet appear immediately.

### Computing a textured rock from EBSD

1. **File → Open EBSD file…** — load a `.ctf` scan.
2. The **EBSD** tab opens automatically showing the phase map.
3. (Optional) Set ROI, apply MAD filter, downsample.
4. In the *VRH Averaging* section, map each detected phase to a material in the database.
5. Set P, T, averaging method (Voigt / Reuss / Hill).
6. Click **Compute VRH** — the textured elastic tensor and *Vp* stereonet appear.
7. Switch to the **Grains** tab → click **Reconstruct grains** for shape and orientation analysis.
8. Switch to **ODF & PF** → **Compute ODF** → plot sections, pole figures, IPF map.

### Computing a modal rock (no EBSD)

1. Open **Modal Rock** tab.
2. Enter phase names and volume fractions (must sum to 100 %).
3. Set P, T.
4. Click **Compute Modal Rock** — Voigt-averaged tensor + stereonet.
5. Switch to the **Hashin-Shtrikman** sub-tab for isotropic upper/lower/VRH bounds.

## Application Layout

```
┌─────────────────────────────────────────────────────────┐
│  Menu: File | Help                                       │
├──────────┬──────────┬──────────┬──────────┬─────────────┤
│ Material │Anisotropy│   EBSD   │Modal Rock│Grains│ODF&PF│
├──────────┴──────────┴──────────┴──────────┴──────────────┤
│                                                           │
│  ┌─────────────────────┐  ┌────────────────────────────┐  │
│  │  Controls / Inputs  │  │    Interactive Plot Area    │  │
│  │   (scrollable)      │  │  (Plotly 2-D or PyVista 3-D)│  │
│  └─────────────────────┘  └────────────────────────────┘  │
├───────────────────────────────────────────────────────────┤
│  Status bar                                               │
└───────────────────────────────────────────────────────────┘
```

Every plot panel has a **Plot Options** group with:

| Control | Description |
|---------|-------------|
| **Colorscale** | Viridis, RdBu_r, Plasma, Jet, … (19 options) |
| **v min / v max** | Clamp the colour range (leave at 0 for auto) |
| **Point size** | Scatter-point radius in pixels |
| **Colorbar** | Toggle the colorbar on/off |

## Tabs — User Guide

### 1. Material Database Tab

**Purpose:** Browse the built-in mineral database and compute elastic properties at arbitrary pressure and temperature.

#### Left panel

| Control | Description |
|---------|-------------|
| Phase table | Searchable list of all phases with crystal system, primary phase flag, and density ρ (g/cm³) |
| Pressure | 0 – 300 GPa (step 0.01 GPa) |
| Temperature | 0 – 3000 K (step 0.1 K) |
| **Compute anisotropy** | Fetch *C*ᵢⱼ at P,T, hand off to AnisotropyBackend |

#### Right panel

| Widget | Description |
|--------|-------------|
| Voigt matrix (6×6) | Single-crystal stiffness in GPa |
| Anisotropy values | Scalar metrics (see table below) |
| Stereonet options | Colorscale, v min/max, point size, colorbar |
| Vp stereonet | Upper-hemisphere stereographic projection of *Vp* (km/s) |

#### Anisotropy scalar metrics

| Metric | Definition |
|--------|-----------|
| **P-wave anisotropy** | `(Vp_max − Vp_min) / Vp_mean × 100` % |
| **S1-wave anisotropy** | `(Vs1_max − Vs1_min) / Vs1_mean × 100` % |
| **S2-wave anisotropy** | `(Vs2_max − Vs2_min) / Vs2_mean × 100` % |
| **AVS** | Max shear-wave splitting `(Vs1 − Vs2) / Vs1_mean × 100` % |
| **Max ΔVs** | Maximum shear-wave splitting in m/s |
| **Vp/Vs1** | Ratio at mean direction |
| **Mean Vp, Vs1, Vs2** | Orientational averages |

#### Pressure–temperature derivatives

Elastic moduli are corrected from their reference values using:

```
Cᵢⱼ(P,T) = Cᵢⱼ⁰ + (dCᵢⱼ/dP)·P + (dCᵢⱼ/dT)·(T − 300 K)
```

where `dCᵢⱼ/dP` and `dCᵢⱼ/dT` are per-component derivatives stored in the database JSON files.

### 2. Anisotropy Tab

**Purpose:** Manually enter a Voigt stiffness matrix and density, compute the full velocity surface, and explore anisotropy interactively in 2-D and 3-D.

#### Left panel

| Control | Description |
|---------|-------------|
| Voigt table (6×6) | Editable cells — enter your *C*ᵢⱼ in GPa |
| Density | 100 – 20 000 kg/m³ |
| n θ | Polar-angle resolution of velocity grid (10 – 300) |
| n φ | Azimuthal-angle resolution (10 – 600) |
| **Compute anisotropy** | Runs in a background thread; progress bar shown |
| Scalar metrics | Vp/Vs ratios, AVS, anisotropy %, Voigt–Reuss–Hill velocities |

#### Right panel

| Control | Description |
|---------|-------------|
| Display combo | Choose *vp* / *vs1* / *vs2* / *avs* / *vpvs1* / *vpvs2* |
| Stereonet options | Colorscale, v min/max, point size, colorbar |
| **2-D stereonet** tab | Upper-hemisphere Plotly scatter coloured by chosen scalar |
| **3-D surface** tab | PyVista velocity surface (sphere scaled by velocity) |

#### Velocity computation

The Christoffel tensor for propagation direction **n** is:

```
Tᵢₖ = Cᵢⱼₖₗ nⱼ nₗ
```

Eigenvalues λ of **T** give phase velocities  `V = √(λ/ρ)`. The three eigenvalues correspond to *Vp* (P-wave), *Vs*₁ and *Vs*₂ (fast and slow S-waves).

#### 3-D surface

The velocity surface is a sphere whose radius in each direction is scaled by the velocity in that direction. Colour maps the same scalar. Rotate, zoom, and pan with the mouse in the PyVista viewport.

### 3. EBSD Tab

**Purpose:** Load, visualise, filter, and process raw EBSD data. Compute the orientation-averaged (textured) elastic tensor via Voigt, Reuss, or Hill averaging.

#### Workflow

```
Load file → View phase map → (Optional) Filter + ROI
    → Map phases to minerals → Set P, T
    → Compute VRH → View stereonet / 3-D surface
```

#### Loading a file

Go to **File → Open EBSD file…** or drag-and-drop. Supported formats:

| Extension | Description |
|-----------|-------------|
| `.ctf` | Oxford HKL Channel 5 text format |
| `.ang` | TSL/EDAX format |
| `.osc` | Oxford binary |
| `.txt` | Generic tab-separated Euler angles |

After loading, the **EBSD summary** panel shows a MTEX-style report with phase percentages, scan dimensions, and step sizes.

#### Phase map

Each phase is assigned a colour from the default palette. Click the coloured button next to each phase name to open a colour picker and choose any custom colour.

Hover over the map to see (x, y, phase) coordinates in the Plotly tooltip.

#### Pre-processing

| Control | Description |
|---------|-------------|
| **MAD filter** | Remove points with Mean Angular Deviation > threshold (default 0.7°) — filters poorly indexed points |
| **Downsample factor** | Keep every *n*-th point (reduces memory for large maps) |

#### Region of Interest (ROI)

Set X min/max and Y min/max spinboxes to clip the dataset to a rectangular subregion. Click **Clip ROI** to apply, **Reset** to restore the full map.

#### VRH Texture Averaging

1. Each detected phase is listed with a **Material** dropdown — select the matching mineral from the database.
2. Set **Pressure** (GPa) and **Temperature** (K).
3. Choose **Averaging method**:
   - **Voigt** — upper bound (strain compatibility assumed); stiffest result
   - **Reuss** — lower bound (stress compatibility assumed); softest result
   - **Hill** — arithmetic mean of Voigt and Reuss; most commonly used
4. Click **Compute VRH** — computation runs in a background thread.

The averaged tensor and density are forwarded to AnisotropyBackend. The *Vp* stereonet and 3-D surface update automatically.

#### Mathematical basis

For *N* EBSD points, the Voigt average is:

```
<Cᵢⱼₖₗ> = (1/N) Σₙ Rᵢₐ(φ₁,Φ,φ₂)ₙ · Rⱼᵦ · Rₖ꜀ · Rₗ꜁ · Cₐᵦ꜀꜁
```

where **R**(φ₁, Φ, φ₂) is the Bunge-convention Euler rotation matrix and *C*ₐᵦ꜀꜁ is the single-crystal tensor in the crystal reference frame.

### 4. Modal Rock Tab

**Purpose:** Compute the effective elastic properties of a multi-mineral rock assemblage without EBSD data (random texture assumed) using either Voigt averaging or Hashin-Shtrikman bounds.

#### Sub-tab 1 — Anisotropic Modal Rock

Assumes **random orientation** (isotropic texture) for each phase and applies a Voigt average of the single-crystal stiffness matrices weighted by volume fraction.

| Control | Description |
|---------|-------------|
| Phase table | Phase name (dropdown) + volume fraction (%) |
| Pressure | GPa |
| Temperature | K |
| **Compute Modal Rock** | Runs Voigt average; shows Voigt matrix, density, stereonet |

> Volume fractions are normalised internally — they do not need to sum exactly to 100 %.

#### Sub-tab 2 — Hashin-Shtrikman Bounds

Computes **isotropic** upper and lower elastic bounds for the composite using the Hashin-Shtrikman variational formulation. Also shows Voigt–Reuss–Hill averages.

| Output | Description |
|--------|-------------|
| **HS Upper** | Stiffest possible isotropic composite (all grains perfectly stiff) |
| **HS Lower** | Softest possible isotropic composite |
| **VRH** | Voigt, Reuss, Hill isotropic averages |
| Bar chart | Grouped bars of Vp / Vs for Lower / VRH / Upper bounds |

**Seismic velocity outputs** per bound:

| Variable | Description |
|----------|-------------|
| *Vp* | P-wave velocity (km/s) |
| *Vs* | S-wave velocity (km/s) |
| *K* | Adiabatic bulk modulus (GPa) |
| *G* | Shear modulus (GPa) |
| ρ | Density (kg/m³) |

### 5. Grains Tab

**Purpose:** Reconstruct individual grains from EBSD data and perform comprehensive shape, orientation, and boundary analysis.

> **Prerequisite:** An EBSD file must be loaded in the EBSD tab first.

#### Sub-tab 1 — Reconstruction

| Control | Description |
|---------|-------------|
| Misorientation threshold | Max allowed intragranular misorientation in degrees (default 10°). Higher → fewer, larger grains |
| Min pixels | Grains with fewer pixels than this are merged into neighbours |
| **Reconstruct grains** | Runs grain reconstruction; progress bar shown |

**Outputs:**
- **Grain ID map** — each pixel coloured by its grain ID (Turbo colorscale)
- **Grain size map** — each pixel coloured by the area of its parent grain

#### Sub-tab 2 — Shape Parameters

After reconstruction, click **Compute shape** to calculate per-grain shape metrics.

| Metric | Description |
|--------|-------------|
| **Area** | Grain area in µm² (or pixels²) |
| **Aspect ratio** | Major / minor ellipse axis |
| **GOS** | Grain Orientation Spread (°) — mean misorientation of all pixels from the grain mean |
| **Ellipse fit** | Semi-axes, orientation angle, centre |
| **Convex hull** | Solidity = grain area / hull area |
| **Feret diameter** | Caliper max/min diameters |

Visualisations:
- Histograms of area, aspect ratio, and GOS distributions
- Scatter: grain size vs. aspect ratio, coloured by GOS

#### Sub-tab 3 — Orientation

| Analysis | Description |
|----------|-------------|
| **KAM** | Kernel Average Misorientation — mean misorientation between each pixel and its neighbours (within `max_angle` threshold). Highlights deformation gradients |
| **GROD** | Grain Reference Orientation Deviation — misorientation of each pixel from its grain mean. Highlights intragranular strain |

Both are rendered as coloured scatter maps. Use the plot-options panel to set the colorscale and range.

#### Sub-tab 4 — Grain Boundaries

| Analysis | Description |
|----------|-------------|
| **Boundary map** | Boundary segments coloured by misorientation angle |
| **Misorientation distribution** | Histogram of boundary misorientation angles. A vertical line marks the LAGB/HAGB threshold |
| **CSL detection** | Identifies Coincidence Site Lattice boundaries (Σ3 = twins, Σ5, Σ7, …) within a tolerance |
| **Twist / Tilt** | Decomposes each boundary segment into twist and tilt components |

The boundary map legend:
- **Grey** — Low-angle grain boundaries (LAGB, < 15° by default)
- **Blue** — High-angle grain boundaries (HAGB, ≥ 15°)
- **Red** — Σ3 twin boundaries (60° / ⟨111⟩ for cubic phases)

#### Sub-tab 5 — Triple Points

Displays the locations of triple-junction points (pixels where three or more grains meet) overlaid on the grain ID map as red × markers.

#### Sub-tab 6 — Export

| Format | Description |
|--------|-------------|
| **CSV** | Grain summary table (ID, phase, area, Euler angles, shape metrics) |
| **Neper ORI** | Orientation file for Neper polycrystal generator |
| **Neper TESS** | Tessellation geometry for Neper |
| **CTF** | Re-exports processed data as an HKL `.ctf` file |

### 6. ODF & PF Tab

**Purpose:** Compute the Orientation Distribution Function (ODF) from EBSD Euler angles and visualise it as Euler sections, pole figures, and inverse pole figure maps.

> **Prerequisite:** An EBSD file must be loaded.

#### ODF Computation

| Control | Description |
|---------|-------------|
| **Phase index** | Which phase to compute the ODF for (1-based) |
| **Crystal symmetry** | All 11 Laue groups supported (see table below) |
| **Sample symmetry** | None or Orthorhombic (orthotropic sample) |
| **Kernel halfwidth** | Smoothing kernel width in degrees (default 10°). Smaller → sharper but noisier |
| **Grid resolution** | ODF evaluation grid spacing in degrees (default 5°) |
| **Max orientations** | Cap on the number of EBSD points used (random subsample if exceeded) |
| **Compute ODF** | Runs kernel-density estimation on SO(3) |

#### Supported Crystal Symmetries

| Laue Group | Crystal System | Schoenflies | Proper ops |
|-----------|---------------|-------------|------------|
| `-1` | Triclinic | Cᵢ | 1 |
| `2/m` | Monoclinic | C₂ₕ | 2 |
| `mmm` | Orthorhombic | D₂ₕ | 4 |
| `4/m` | Tetragonal (low) | C₄ₕ | 4 |
| `4/mmm` | Tetragonal | D₄ₕ | 8 |
| `-3` | Trigonal (low) | S₆ | 3 |
| `-3m` | Trigonal / Rhombohedral | D₃d | 6 |
| `6/m` | Hexagonal (low) | C₆ₕ | 6 |
| `6/mmm` | Hexagonal | D₆ₕ | 12 |
| `m-3` | Cubic (low) | Tₕ | 12 |
| `m-3m` | Cubic | Oₕ | 24 |

Aliases are accepted: `olivine` → `mmm`, `cubic` → `m-3m`, `Ci` → `-1`, etc.

#### ODF Sections

Shows the ODF as a grid of 2-D contour plots at fixed φ₂ (or φ₁ or σ) values.

| Control | Description |
|---------|-------------|
| Section type | `phi2 sections`, `phi1 sections`, `sigma sections` |
| Section values | Comma-separated list of angles in degrees (e.g. `0, 15, 30, 45, 60, 75, 90`) |
| Section plot options | Colorscale, v min/max (m.r.d.), colorbar |

Values are in **multiples of random distribution (m.r.d.)**. A value of 1 indicates a random texture. Values > 1 indicate preferred orientations.

#### Pole Figure

Renders a stereographic projection of a set of crystallographic directions specified by Miller indices (h k l).

| Control | Description |
|---------|-------------|
| h k l | Miller indices of the reflection (e.g. 1 0 0, 1 1 1) |
| KDE bandwidth | Kernel smoothing bandwidth in degrees (default 7.5°) |
| PF plot options | Colorscale, v min/max, point size, colorbar |

The plot shows:
- **Contour fill** — KDE density in m.r.d.
- **Raw scatter** — individual EBSD pole positions (semi-transparent)
- **N / S / E / W** reference labels

#### Inverse Pole Figure (IPF) Map

Colours each EBSD point by its crystallographic direction aligned with a chosen sample direction:

| Sample direction | Meaning |
|-----------------|---------|
| Z (ND) | Normal direction (out of plane) |
| X (RD) | Rolling direction |
| Y (TD) | Transverse direction |

The right sub-panel shows the IPF colour triangle (fundamental domain of the crystal symmetry group).

#### Texture Components

Computes the volume fraction of standard texture components within a user-defined angular tolerance (default 15°). Components are shown as a horizontal bar chart with descriptions.

#### ODF Properties

Displays quantitative texture descriptors:

| Property | Description |
|----------|-------------|
| **Texture index J** | `∫ f(g)² dg` — J = 1 for random, J > 1 for textured |
| **Entropy S** | `−∫ f(g) ln f(g) dg` — S = 0 for random, S < 0 for textured |
| **ODF max** | Peak value in m.r.d. |

#### Export Formats

| Format | Use |
|--------|-----|
| **POPLA `.epf`** | EPSC (Elasto-Plastic Self-Consistent) codes |
| **MTEX-compatible text** | Re-import into MTEX for further analysis |
| **VPSC orientations** | Los Alamos VPSC crystal plasticity code |

## Python API Reference

All core functionality is available as a Python library independent of the GUI. Import from the `santex` package:

```python
import santex
from santex.material.material import Material
from santex.ebsd.ebsd import EBSD
from santex.grains import Grain2D, calc_grains
from santex.anisotropy.anisotropy import Anisotropy
from santex.odf.odf import ODF
from santex.isotropy.isotropy import Isotropy
```

### `material/material.py`

```python
class Material(database_path, database_path2, database_path3)
```

Loads three JSON databases (ambient constants, pressure derivatives, temperature derivatives) and provides interpolated stiffness matrices at arbitrary P, T.

#### Methods

```python
get_voigt_matrix(phase: str) -> np.ndarray  # shape (6, 6), GPa
```
Return the ambient (P=0, T=300 K) Voigt stiffness matrix.

```python
voigt_high_PT(phase: str, pressure: float, temp: float) -> np.ndarray
```
Return P,T-corrected Voigt matrix using linear derivatives.
- `pressure` in GPa
- `temp` in K

```python
load_density(phase: str, pressure=None, temperature=None) -> float
```
Return density in kg/m³ at the given conditions.

```python
modal_rock(
    minerals: list[str],
    fractions: list[float],
    pressure: float,
    temp: float
) -> tuple[np.ndarray, float]
```
Voigt average of multiple minerals.  Returns `(cij_GPa, rho_kg_m3)`.

**Example:**
```python
from santex.material.material import Material

mat = Material(db_path, deriv_P_path, deriv_T_path)
cij = mat.voigt_high_PT("forsterite", pressure=3.0, temp=1000.0)
rho = mat.load_density("forsterite", pressure=3.0, temperature=1000.0)
```

### `ebsd/ebsd.py`

```python
class EBSD(filename: str)
```

Load and process EBSD data from `.ctf`, `.ang`, `.osc`, or `.txt` files.

#### Methods

```python
phases() -> list[tuple[int, str, float]]
```
Returns `[(phase_id, name, percentage), ...]`.

```python
get_euler_angles(phase: int, data=None) -> pd.DataFrame
```
Returns a DataFrame with columns `['Euler1', 'Euler2', 'Euler3']` (degrees) for the requested phase.

```python
get_voigt_reuss_hill(
    cij: np.ndarray,          # (6,6) GPa
    euler_angles: pd.DataFrame,
    density: float,           # kg/m³
    method: str = "voigt"     # "voigt" | "reuss" | "hill"
) -> tuple[np.ndarray, float]
```
Returns the orientation-averaged `(cij_GPa, rho_kg_m3)`.

```python
filter_MAD(data: pd.DataFrame, threshold: float = 0.7) -> pd.DataFrame
```
Remove rows where MAD > threshold.

```python
downsample_EBSD(factor: int, data: pd.DataFrame) -> pd.DataFrame
```

**Example:**
```python
from santex.ebsd.ebsd import EBSD
from santex.material.material import Material

ebsd = EBSD("my_scan.ctf")
euler = ebsd.get_euler_angles(phase=1)
euler = ebsd.filter_MAD(euler)

mat = Material(...)
cij0 = mat.voigt_high_PT("forsterite", 0, 300)
rho  = mat.load_density("forsterite")

cij_avg, rho_avg = ebsd.get_voigt_reuss_hill(cij0, euler, rho, method="hill")
```

### `grains/grains.py`

```python
def calc_grains(
    data: pd.DataFrame,
    threshold_deg: float = 10.0,
    min_pixels: int = 2
) -> Grain2D
```

Reconstruct grains from EBSD data using a misorientation flood-fill algorithm.

```python
class Grain2D(data, grain_ids, all_pairs=None)
```

Container for a set of reconstructed grains.

#### Key properties (all cached, per-grain arrays)

| Property | Shape | Description |
|----------|-------|-------------|
| `grain_id_list` | (N_grains,) | Integer IDs |
| `area` | (N_grains,) | µm² |
| `aspect_ratio` | (N_grains,) | Major/minor |
| `diameter` | (N_grains,) | Equivalent circular diameter |
| `mean_euler` | (N_grains, 3) | Mean Euler angles (deg) |
| `gos` | (N_grains,) | Grain orientation spread (deg) |

#### Key methods

```python
grain2d.fit_ellipse()        -> dict   # semi-axes, angle, centre
grain2d.feret_diameter()     -> dict   # max, min, mean
grain2d.convex_hull_props()  -> dict   # area, perimeter, solidity
grain2d.boundary             -> GrainBoundary
grain2d.triple_points        -> TriplePoints
grain2d.summary_df()         -> pd.DataFrame
```

```python
def calc_kam(
    data: pd.DataFrame,
    grain_ids: np.ndarray,
    max_angle_deg: float = 5.0
) -> np.ndarray
```
Returns per-pixel KAM values (degrees).

```python
def calc_grod(
    data: pd.DataFrame,
    grain_ids: np.ndarray
) -> np.ndarray
```
Returns per-pixel GROD values (degrees).

**Example:**
```python
from santex.grains import calc_grains, calc_kam

grains = calc_grains(ebsd_data, threshold_deg=10, min_pixels=3)
print(f"Reconstructed {grains.n_grains} grains")
print(grains.summary_df().head())

kam = calc_kam(ebsd_data, grains.get_grain_ids())
```

### `anisotropy/anisotropy.py`

```python
class Anisotropy(stiffness_matrix: np.ndarray, density: float)
```

Compute phase velocities and anisotropy metrics from a Voigt stiffness matrix.
- `stiffness_matrix` — shape (6, 6), **Pa** (not GPa)
- `density` — kg/m³

#### Methods

```python
phase_velocity() -> tuple[np.ndarray, np.ndarray, np.ndarray]
```
Returns `(Vp, Vs1, Vs2)` over a spherical grid, in m/s.

```python
anisotropy_values() -> dict
```
Returns a dict of scalar metrics (see Material Database tab for list).

```python
voigt_velocity() -> tuple   # (Vp, Vs1, Vs2, Vbulk) in m/s
reuss_velocity() -> tuple
hill_velocity()  -> tuple
```

**Example:**
```python
from santex.anisotropy.anisotropy import Anisotropy
import numpy as np

# Forsterite stiffness (GPa → Pa)
cij = np.array([
    [320.5, 68.2, 71.6,  0,    0,    0  ],
    [ 68.2,195.5, 76.8,  0,    0,    0  ],
    [ 71.6, 76.8,233.5,  0,    0,    0  ],
    [  0,    0,    0,   64.6,  0,    0  ],
    [  0,    0,    0,    0,   77.0,  0  ],
    [  0,    0,    0,    0,    0,   78.7],
]) * 1e9

anis = Anisotropy(cij, density=3355.0)
metrics = anis.anisotropy_values()
print(f"P-wave anisotropy: {metrics['P wave anisotropy percent']:.2f} %")
print(f"AVS max: {metrics['Max vs anisotropy percent']:.2f} %")
```

### `odf/odf.py`

```python
class ODF(
    euler_deg: np.ndarray,    # (N, 3) Euler angles in degrees
    crystal_symmetry: str,    # any alias accepted (e.g. "mmm", "m-3m", "olivine")
    halfwidth_deg: float = 10.0
)
```

Estimate the ODF by kernel-density estimation on SO(3) using the de la Vallée-Poussin kernel.

#### Methods

```python
phi2_sections(
    phi2_values_deg: list[float],
    resolution_deg: float = 5.0
) -> dict[float, tuple[np.ndarray, np.ndarray, np.ndarray]]
```
Returns `{phi2_angle: (phi1_grid, PHI_grid, odf_values)}`.

```python
get_pole_figure_kde(
    hkl: tuple[int, int, int],
    sample_symmetry: str,
    bandwidth_deg: float = 7.5
) -> tuple[np.ndarray, np.ndarray, np.ndarray]
```
Returns `(xy_grid, density, xy_scatter)` for stereographic projection.

```python
texture_index(resolution_deg: float = 10.0) -> float
entropy(resolution_deg: float = 10.0) -> float
odf_max(resolution_deg: float = 5.0) -> float
```

```python
component_analysis(max_angle_deg: float = 15.0) -> list[dict]
```
Returns list of `{"name": str, "fraction": float, "description": str}`.

```python
random_sampling(n: int = 1000) -> np.ndarray  # (n, 3) Euler angles, degrees
```

**Example:**
```python
from santex.odf.odf import ODF
import numpy as np

euler = ebsd.get_euler_angles(phase=1).values  # (N, 3)
odf = ODF(euler, crystal_symmetry="mmm", halfwidth_deg=10)

J = odf.texture_index()
S = odf.entropy()
print(f"Texture index J = {J:.3f}  (1 = random)")
print(f"Entropy S = {S:.3f}  (0 = random)")

sections = odf.phi2_sections([0, 15, 30, 45, 60, 75, 90], resolution_deg=5)
```

### `isotropy/isotropy.py`

```python
class Isotropy()
```

Compute isotropic seismic velocities using the Birch-Murnaghan equation of state (adapted from Hacker & Abers, 2004).

#### Methods

```python
get_available_phases() -> list[str]
```

```python
calculate_seismic_properties(
    phase: str,
    temperature: float,     # °C
    pressure: float,        # GPa
    ref_density=None,
    return_vp_vs_vbulk=False,
    return_aktout=False
) -> tuple
```

Returns `(density, K_adiabatic, G_shear)` or `(density, K, G, Vp, Vs, Vbulk)`.

**Example:**
```python
from santex.isotropy.isotropy import Isotropy

iso = Isotropy()
rho, K, G, Vp, Vs, Vb = iso.calculate_seismic_properties(
    "forsterite", temperature=800, pressure=3.0,
    return_vp_vs_vbulk=True
)
print(f"Vp = {Vp/1000:.3f} km/s  Vs = {Vs/1000:.3f} km/s")
```

### `tensor/tensor.py`

```python
class Tensor(data=None)
```

Utilities for Voigt notation and Euler-angle rotations.

#### Methods

```python
tensor_to_voigt(tensor: np.ndarray) -> np.ndarray  # (3,3,3,3) → (6,6)
voigt_to_tensor(voigt: np.ndarray) -> np.ndarray   # (6,6) → (3,3,3,3)
euler_to_rotation(phi1, phi, phi2) -> np.ndarray   # degrees → (3,3) rotation matrix
rotate_tensor(tensor, phi1, phi, phi2) -> np.ndarray  # (3,3,3,3) → rotated
```

The Euler rotation uses Bunge (ZXZ) convention:

```
R = Rz(φ₁) · Rx(Φ) · Rz(φ₂)
```

### `grains/symmetry.py`

Provides all 32 crystallographic point groups organised into the 11 Laue groups, plus MTEX-compatible aliases.

```python
from santex.grains.symmetry import (
    get_symmetry,       # get_symmetry("mmm") -> (N, 3, 3) rotation matrices
    get_laue_group,     # get_laue_group("olivine") -> "mmm"
    get_info,           # get_info("m-3m") -> full metadata dict
    list_symmetries,    # list_symmetries() -> all 107 aliases
    list_laue_groups,   # list_laue_groups() -> 11 Laue-group dicts
    SYMMETRY,           # dict mapping every alias to rotation matrices
    CRYSTAL_SYMMETRY_INFO,  # full metadata per Laue group
)
```

#### Accepted aliases (examples)

| Input | Resolved Laue group |
|-------|-------------------|
| `"olivine"`, `"forsterite"`, `"enstatite"` | `mmm` |
| `"cubic"`, `"m-3m"`, `"Oh"`, `"43"` | `m-3m` |
| `"hexagonal"`, `"6/mmm"`, `"D6h"` | `6/mmm` |
| `"trigonal"`, `"-3m"`, `"D3d"` | `-3m` |
| `"tetragonal"`, `"4/mmm"`, `"D4h"` | `4/mmm` |
| `"monoclinic"`, `"2/m"`, `"C2h"` | `2/m` |
| `"triclinic"`, `"-1"`, `"Ci"` | `-1` |

## Architecture

```
santex_app/
├── main.py                        ← QApplication entry point
└── frontend/
    ├── main_window.py             ← MainWindow (QMainWindow)
    ├── tabs/
    │   ├── material_tab.py        ← MaterialTab
    │   ├── anisotropy_tab.py      ← AnisotropyTab
    │   ├── ebsd_tab.py            ← EBSDTab + _PlotOptions
    │   ├── grains_tab.py          ← GrainsTab
    │   ├── modal_rock_tab.py      ← ModalRockTab
    │   └── odf_tab.py             ← ODFTab
    └── widgets/
        ├── plotly_widget.py       ← PlotlyWidget (QWebEngineView)
        └── pyvista_widget.py      ← PyVistaWidget (BackgroundPlotter)

santex_app/backend/
├── material_backend.py            ← MaterialBackend
├── ebsd_backend.py                ← EBSDBackend
├── grains_backend.py              ← GrainsBackend
├── anisotropy_backend.py          ← AnisotropyBackend
└── odf_backend.py                 ← ODFBackend

santex/                            ← Core library (GUI-independent)
├── material/material.py
├── ebsd/ebsd.py
├── grains/
│   ├── grain2d.py
│   ├── symmetry.py
│   └── ...
├── anisotropy/anisotropy.py
├── odf/odf.py
├── isotropy/isotropy.py
└── tensor/tensor.py
```

### Design principles

**Layered architecture**

```
GUI Widgets  →  Tab logic  →  Backends  →  Core library
(Qt/Plotly)     (PyQt5)        (thin)       (pure Python/numpy)
```

- **Core library** (`santex/`) has zero Qt dependency. Can be used from notebooks or scripts.
- **Backends** (`backend/`) are thin wrappers that translate between the core library's data types and what the GUI needs. No plotting, no Qt widgets.
- **Tabs** (`frontend/tabs/`) own the UI layout and call backends. Long computations run in `QThread` workers to keep the GUI responsive.
- **Widgets** (`frontend/widgets/`) are reusable rendering components.

### Worker thread pattern

Every long computation (grain reconstruction, ODF, VRH, velocity surface) follows the same pattern:

```python
class _Worker(QThread):
    finished = pyqtSignal(object)
    error    = pyqtSignal(str)

    def run(self):
        try:
            result = backend.compute_something(...)
            self.finished.emit(result)
        except Exception as e:
            self.error.emit(str(e))

# In the tab:
worker = _Worker(backend, ...)
worker.finished.connect(self._on_done)
worker.error.connect(self._on_error)
worker.start()
```

### PlotlyWidget rendering pipeline

```
go.Figure object
    → pio.to_html(include_plotlyjs=False)
    → inject CSS-insertRule shim (fixes Qt 5.15 Chromium compatibility)
    → inject <script src="file:///TEMPDIR/plotly.min.js">
    → write to TEMPDIR/plot_<id>.html
    → QWebEngineView.load(QUrl.fromLocalFile(...))
```

The `plotly.min.js` is copied into the same temp directory as the HTML at startup so both share the same `file://` origin, avoiding Qt WebEngine's cross-origin restrictions.

## File Formats

### Input formats

| Format | Extension | Notes |
|--------|-----------|-------|
| HKL Channel 5 | `.ctf` | Text, tab-delimited Euler angles + phase |
| TSL/EDAX | `.ang` | Text, Euler angles + confidence index |
| Oxford binary | `.osc` | Binary format |
| Generic text | `.txt` | Tab-separated with header |

#### `.ctf` file structure (example header)

```
Channel Text File
Prj   my_project
Author  John Doe
JobMode Grid
XCells  400
YCells  300
XStep   0.1
YStep   0.1
AcqE1   0
AcqE2   0
AcqE3   0
Euler angles refer to Sample Coordinate system (CS0)!
Phases 2
3.956;3.956;3.956 90;90;90 Spinel ... 
4.756;10.227;5.994 90;105.1;90 Monoclinic ...
Phase X Y Bands Error Euler1 Euler2 Euler3 MAD BC BS
1 0.1 0.1 9 0 125.3 58.2 312.1 0.41 188 255
```

### Output formats

| Format | Extension | Producer | Consumer |
|--------|-----------|---------|---------|
| POPLA | `.epf` | ODF export | EPSC, popLA |
| MTEX text | `.txt` | ODF export | MTEX toolbox |
| VPSC orientations | `.txt` | ODF export | Los Alamos VPSC |
| Grain summary | `.csv` | Grains export | Excel, Python, R |
| Neper ORI | `.ori` | Grains export | Neper polycrystal |
| Neper TESS | `.tess` | Grains export | Neper FEM |
| HKL CTF | `.ctf` | EBSD / Grains export | AZtec, MTEX, HKL |

## Keyboard Shortcuts

| Shortcut | Action |
|----------|--------|
| `Ctrl+O` | Open EBSD file |
| `Ctrl+Q` | Quit |

### Plotly chart interactions (mouse)

| Action | Effect |
|--------|--------|
| Click + drag | Pan |
| Scroll wheel | Zoom |
| Double-click | Reset zoom |
| Hover | Show tooltip (X, Y, value) |
| Box select | Zoom to region |
| Camera icon (toolbar) | Export to SVG |

### PyVista viewport (mouse)

| Action | Effect |
|--------|--------|
| Left drag | Rotate |
| Right drag / scroll | Zoom |
| Middle drag | Pan |
| `r` key | Reset camera |

## Troubleshooting

### App does not start — `ModuleNotFoundError: PyQt5.QtWebEngineWidgets`

The pip and conda PyQt5 packages use different Qt DLL naming conventions on Windows. Replace the conda PyQt5 with the pip version:

```bash
conda activate santex
pip install --force-reinstall PyQt5 PyQtWebEngine
```

### Plots show blank / "Plotly is not defined"

Occurs when Qt WebEngine cannot load `plotly.min.js` due to cross-origin restrictions. Fixed in v1.2.3+ — the bundled JS is now copied into the same temp directory as the HTML. If you still see this:

1. Check that `plotly` is installed: `python -c "import plotly; print(plotly.__version__)"`
2. Verify the temp dir is writable: `python -c "import tempfile; print(tempfile.mkdtemp())"`

### Grain reconstruction is very slow

- Reduce the EBSD dataset first: apply **MAD filter** and use a **downsample factor** of 5–10.
- Lower **n θ** and **n φ** in the Anisotropy tab.
- Increase **Min pixels** in the Grains tab to ignore tiny grains.

### VTK / PyVista 3-D surface does not render

The 3-D surface requires a working OpenGL context. In headless / remote-desktop environments this may fail. The error `wglMakeCurrent failed` is harmless if it appears only in the terminal — the 2-D stereonet in the Anisotropy tab is a complete alternative.

### ODF computation gives wrong result

- Check that the correct **crystal symmetry** is selected. Olivine → `mmm`, quartz → `-3m`, calcite → `-3m`, feldspar → `-1`.
- A **halfwidth** that is too large (> 15°) will over-smooth and flatten the ODF.
- Ensure the **phase index** matches the phase you want (check the EBSD summary for phase IDs).

### EBSD file loads but shows no data

- Confirm the file is a valid CTF / ANG / OSC — open it in a text editor and check for the expected header.
- For `.ctf` files from older HKL software, the column order may differ; check `X`, `Y`, `Phase`, `Euler1/2/3`, `MAD` columns are present.

### Pressure–temperature derivatives seem wrong

The database covers:
- Pressure: 0 – ~300 GPa (extrapolation beyond measured range is flagged)
- Temperature: 0 – ~3000 K

Derivatives are first-order linear. For large P,T excursions, results are approximate.

## References

- Singh *et al.* (2025). *SAnTex: Seismic Anisotropy from Texture*. Journal of Open Source Software.
- Hielscher & Schaeben (2008). *A novel pole figure inversion method: specification of the MTEX algorithm*. J. Appl. Cryst. 41, 1024–1037.
- Hacker & Abers (2004). *Subduction Factory 3: An Excel worksheet and macro for calculating the densities, seismic wave speeds, and H₂O contents of minerals and rocks at pressure and temperature*. Geochem. Geophys. Geosyst. 5.
- Mainprice (1990). *A FORTRAN program to calculate seismic anisotropy from the lattice preferred orientation of minerals*. Computers & Geosciences 16(3), 385–393.
- Bunge (1982). *Texture Analysis in Materials Science*. Butterworths, London.
- International Tables for Crystallography, Volume A (Space-Group Symmetry). IUCr / Springer.

*Documentation generated for SAnTex v1.2.3 — last updated May 2026.*
