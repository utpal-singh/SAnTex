# Introduction

**SAnTex** (Seismic Anisotropy from Texture) is an open-source Python application and library for computing the full elastic tensor of rocks from crystallographic texture data and calculating the resulting seismic wave velocities and anisotropy.

It combines a material database of single-crystal elastic constants, EBSD data processing, grain reconstruction, Orientation Distribution Function (ODF) analysis, and interactive 2-D / 3-D visualisation — all in a single desktop GUI.

## What SAnTex does

| Step | Input | Output |
|------|-------|--------|
| **1. Single-crystal database** | Phase name, P, T | Voigt stiffness matrix *C*ᵢⱼ (GPa) |
| **2. Texture averaging** | EBSD Euler angles + *C*ᵢⱼ | Polycrystal tensor (Voigt / Reuss / Hill) |
| **3. Velocity computation** | Polycrystal *C*ᵢⱼ, density ρ | *Vp*, *Vs*₁, *Vs*₂ surfaces (km/s) |
| **4. Anisotropy metrics** | Velocity surfaces | AVS, P-wave anisotropy, Vp/Vs ratios |
| **5. Modal rock** | Volume fractions + P, T | Rock-average tensor + isotropic HS bounds |

## Key features

- **Material database** — 55+ rock-forming minerals with pressure and temperature derivatives (first-order linear EOS)
- **EBSD processing** — load `.ctf` / `.ang` / `.osc` files, filter by MAD, define spatial ROI, downsample
- **Grain reconstruction** — misorientation-based segmentation, shape analysis (ellipse, Feret, convex hull), KAM, GROD, CSL boundary detection
- **ODF / Pole figures** — kernel-density ODF on SO(3), all 11 Laue-group symmetries, KDE pole figures, IPF maps
- **Interactive plots** — Plotly (2-D stereonets, contour ODF sections, bar charts) + PyVista (3-D velocity surface)
- **Export** — POPLA `.epf`, MTEX-compatible text, VPSC orientation files, grain CSV, Neper TESS/ORI, `.ctf`

## Scientific background

The elastic stiffness tensor **C**ᵢⱼₖₗ relates stress σᵢⱼ to strain εₖₗ via Hooke's Law:

$$\sigma_{ij} = C_{ijkl} \, \epsilon_{kl}$$

Its pressure and temperature dependence is linearised as:

$$C_{ijkl}(p,T) = C_{ijkl}(0,0)
  + \left.\frac{\partial C_{ijkl}}{\partial p}\right|_0 p
  + \left.\frac{\partial C_{ijkl}}{\partial T}\right|_0 T
  + \mathcal{O}(p^2, T^2)$$

For a polycrystal with *N* grains, the Voigt average is:

$$\bar{C}^{\text{Voigt}}_{ijkl} = \frac{1}{N}\sum_{n=1}^{N} R_{ia}R_{jb}R_{kc}R_{ld}\,C_{abcd}$$

where **R** is the rotation matrix for each grain orientation (Bunge ZXZ Euler angles). The Reuss average uses the compliance tensor **S** = **C**⁻¹, and Hill takes the arithmetic mean of the two.

## Architecture overview

```
EBSD file (.ctf / .ang)
    ↓  EBSDBackend
Euler angles per grain
    ↓  AnisotropyBackend.vrh_average()
Polycrystal Cᵢⱼₖₗ (Pa)
    ↓  Christoffel equation
Vp, Vs1, Vs2 surfaces
    ↓  PlotlyWidget / PyVistaWidget
2-D stereonet / 3-D sphere
```

The core library (`santex/`) is completely independent of Qt and can be used from Jupyter notebooks or scripts. The GUI adds backends (thin wrappers) and tabs (Qt widgets).

::: tip Next step
Head to [Installation](/guide/installation) to set up the environment.
:::
