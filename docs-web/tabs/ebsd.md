# EBSD Tab

Load, filter, and visualise EBSD data, then compute the Voigt–Reuss–Hill averaged elastic tensor for a textured polycrystal.

## Panel layout

| Panel | Contents |
|-------|----------|
| **Left top** | File path display + Open button |
| **Left** | ROI spinboxes (X/Y start/end) |
| **Left** | MAD threshold, downsample factor |
| **Left** | Phase-to-mineral mapping table |
| **Left** | P, T, averaging method (V/R/H) |
| **Left** | *Compute VRH* button |
| **Right** | Plot type selector + colour options |
| **Right** | Phase-colour pickers |
| **Right** | Interactive Plotly map |

## Supported file formats

| Format | Extension | Source |
|--------|-----------|--------|
| HKL Channel 5 | `.ctf` | Oxford Instruments / Aztec |
| TSL / EDAX | `.ang` | EDAX OIM |
| Oxford binary | `.osc` | Oxford Instruments |
| Generic text | `.txt` | Tab-separated with header |

## Controls reference

| Control | Default | Description |
|---------|---------|-------------|
| ROI (X start / X end) | full | Restrict scan columns |
| ROI (Y start / Y end) | full | Restrict scan rows |
| MAD threshold | 1.0 | Remove pixels with Mean Angular Deviation > threshold |
| Downsample factor | 1 | Keep every n-th pixel (speeds up grain recon + ODF) |
| Phase mapping | — | Map integer phase ID → mineral name |
| Pressure | 0 GPa | For elastic constant lookup |
| Temperature | 300 K | For elastic constant lookup |
| Averaging | Hill | Voigt / Reuss / Hill |
| Plot type | Phase map | Phase map · IPF · Euler · Band Contrast |

## VRH averaging

For a scan with *N* pixels, the Voigt average elastic tensor is:

$$\bar{C}^{\text{Voigt}}_{ijkl} = \frac{1}{N}\sum_{n=1}^{N} R_{ia}^{(n)} R_{jb}^{(n)} R_{kc}^{(n)} R_{ld}^{(n)}\,C_{abcd}$$

where $\mathbf{R}^{(n)}$ is the rotation matrix for pixel *n* (from its Bunge Euler angles $\phi_1, \Phi, \phi_2$).

Reuss uses the compliance **S** = **C**⁻¹:

$$\bar{S}^{\text{Reuss}}_{ijkl} = \frac{1}{N}\sum_{n=1}^{N} R_{ia}^{(n)} R_{jb}^{(n)} R_{kc}^{(n)} R_{ld}^{(n)}\,S_{abcd}$$

Hill (VRH) is the arithmetic mean: $\bar{C}^{\text{Hill}} = \tfrac{1}{2}(\bar{C}^{\text{Voigt}} + (\bar{S}^{\text{Reuss}})^{-1})$.

## Workflow

```
1.  File → Open EBSD file (Ctrl+O)
2.  [Optional] Set ROI  →  reduces dataset
3.  [Optional] MAD ≤ 0.8  →  removes poorly indexed pixels
4.  [Optional] Downsample 5–10×  →  speeds up reconstruction
5.  Map phase IDs to minerals (Phase 1 → "Forsterite", etc.)
6.  Set P, T, averaging method
7.  Click "Compute VRH"
8.  Vp stereonet + scalar metrics appear
```

## Plot types

| Type | Description |
|------|-------------|
| **Phase map** | Each phase coloured by its assigned colour |
| **IPF** | Inverse Pole Figure coloured by crystallographic orientation |
| **Euler** | RGB encoding of φ₁, Φ, φ₂ angles |
| **Band contrast** | Grey-scale EBSD quality map |

## Export

Right-click the phase table header or use **File → Export** for:
- Filtered EBSD data → `.ctf`
- Phase summary → `.csv`
