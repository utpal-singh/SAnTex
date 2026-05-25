# Grains Tab

Reconstruct grains from EBSD data and analyse per-grain shape, orientation, and misorientation statistics.

::: tip Prerequisite
Load an EBSD file in the **EBSD** tab first.
:::

## Panel layout

| Panel | Contents |
|-------|----------|
| **Left** | Misorientation threshold, min pixels |
| **Left** | *Reconstruct grains* button |
| **Left** | Export buttons (CSV, CTF, Neper) |
| **Right** | Plot type selector |
| **Right** | Interactive grain map (Plotly) |
| **Right** | Grain statistics table |

## Controls

| Control | Default | Description |
|---------|---------|-------------|
| Misorientation threshold | 10 ° | Grain boundary detection threshold |
| Min pixels | 5 | Remove grains smaller than this many pixels |
| Plot type | Grain map | See plot types below |

## Grain reconstruction algorithm

1. For every pair of adjacent pixels, compute the misorientation angle.
2. If misorientation ≥ threshold → grain boundary.
3. Connected-component labelling assigns grain IDs.
4. Grains with fewer than *min pixels* are dissolved into neighbours.

## Computed grain properties

| Property | Description |
|----------|-------------|
| Grain ID | Integer label |
| Phase | Phase index |
| N pixels | Size in pixels |
| Area (µm²) | Physical area from step size |
| Centroid X, Y | Centre of mass |
| Mean φ₁, Φ, φ₂ | Mean Euler angles |
| Misorientation | Mean orientation from grain average |
| **KAM** | Kernel Average Misorientation (local plasticity indicator) |
| **GROD** | Grain Reference Orientation Deviation |
| Aspect ratio | Major / minor ellipse axis |
| Ellipse angle | Orientation of major axis |
| Feret max/min | Maximum and minimum caliper diameters |
| CSL type | Coincidence Site Lattice boundary type (Σ3, Σ5, …) |

## Plot types

| Type | Description |
|------|-------------|
| **Grain map** | Each grain a unique colour |
| **KAM** | Kernel Average Misorientation (colour = local strain) |
| **GROD** | Grain Reference Orientation Deviation |
| **Phase** | Phase-coloured grain map |
| **Grain size** | Colour by area |
| **Aspect ratio** | Colour by elongation |

## Export formats

| Format | Button | Contents |
|--------|--------|----------|
| CSV | Export CSV | All grain properties |
| CTF | Export CTF | HKL Channel 5 — grain-average orientations |
| Neper ORI | Export Neper | Euler angles for Neper polycrystal FEM |
| Neper TESS | Export Neper | Voronoi tessellation for Neper |
