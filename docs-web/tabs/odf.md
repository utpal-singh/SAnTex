# ODF & Pole Figures Tab

Compute the Orientation Distribution Function (ODF) from EBSD data, visualise φ₂ sections, plot pole figures, and generate inverse pole figure (IPF) maps.

::: tip Prerequisite
Load and (optionally) filter an EBSD file in the **EBSD** tab first.
:::

## Panel layout

| Panel | Contents |
|-------|----------|
| **Left** | Phase selector, crystal symmetry dropdown |
| **Left** | Halfwidth (°), resolution (°) |
| **Left** | *Compute ODF* button |
| **Left** | Texture index *J*, entropy |
| **Left** | Export buttons (POPLA, MTEX, VPSC) |
| **Right** | Plot type tabs: φ₂ sections · Pole figure · IPF map |
| **Right** | Per-plot options: colorscale, φ₂ angle, pole direction |

## Controls

| Control | Default | Description |
|---------|---------|-------------|
| Phase | 1 | Phase index from EBSD data |
| Crystal symmetry | `mmm` | Laue class (see table below) |
| Halfwidth | 5 ° | KDE bandwidth — smaller = sharper ODF |
| Resolution | 5 ° | ODF grid step size |
| φ₂ section | 0 ° | Which φ₂ = const section to display |
| Pole direction | [100] | Miller indices for pole figure |

## ODF computation

SAnTex uses a **kernel density estimate** (KDE) on SO(3) — the space of 3-D rotations — via the `orix` library:

$$f(g) = \frac{1}{N}\sum_{n=1}^{N} \psi_h(g \cdot g_n^{-1})$$

where $g_n$ are the measured orientations, $\psi_h$ is a de la Vallée Poussin kernel with halfwidth $h$, and $g$ is a point on a uniform SO(3) grid.

## Crystal symmetry — Laue groups

| Laue class | Minerals |
|-----------|---------|
| `-1` | Albite, triclinic feldspars, kaolinite |
| `2/m` | Biotite, muscovite, diopside, hornblende |
| `mmm` | Olivine (Forsterite), enstatite, orthopyroxene |
| `4/m` | — |
| `4/mmm` | Rutile |
| `-3` | — |
| `-3m` | Quartz, calcite |
| `6/m` | — |
| `6/mmm` | — |
| `m-3` | — |
| `m-3m` | Garnet (cubic) |

## ODF properties

| Property | Symbol | Definition |
|----------|--------|------------|
| **Texture index** | J | $J = \int_{SO(3)} f(g)^2 \, dg$ — 1 = random, ∞ = single crystal |
| **Entropy** | S | $S = -\int f(g)\ln f(g)\,dg$ — higher = more random |
| **Max ODF** | — | Peak ODF value |
| **Pfeiffer index** | — | Normalised spread measure |

## Plot types

### φ₂ sections
Constant-φ₂ slices through the 3-D Euler space (φ₁, Φ axes shown). Standard output for orthorhombic samples (olivine, quartz). Multiple φ₂ values plotted in a grid.

### Pole figures
Stereographic projection of a specific crystallographic direction. Upper hemisphere shown by default. Contoured using KDE.

### IPF map
Inverse Pole Figure colour-coded EBSD map: each pixel is coloured by the crystallographic direction parallel to the sample *Z* axis. The colour key (standard crystallographic triangle) is shown.

## Export formats

| Format | File | Consumer |
|--------|------|---------|
| POPLA | `.epf` | EPSC, popLA codes |
| MTEX | `.txt` | MTEX Matlab toolbox |
| VPSC | `.txt` | Los Alamos VPSC |
