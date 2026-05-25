# `santex.grains` — Grain Analysis

Grain reconstruction and shape / misorientation statistics.

## Usage via EBSD

Most grain operations are accessed through the `EBSD` class:

```python
from santex.ebsd import EBSD
ebsd = EBSD("scan.ctf")
df = ebsd.get_ebsd_data()
df_clean = ebsd.filter_MAD(df, 0.8)

phases = ebsd.phases_names()['phase'].tolist()
phases.insert(0, "na")   # phase index 0 = unindexed

df_grains = ebsd.calc_grains(
    df=df_clean,
    threshold=10,       # misorientation threshold (degrees)
    phase_names=phases,
    downsampling_factor=5,
)
```

## `santex.grains.grain2d` — Shape descriptors

Computed automatically during reconstruction and appended as columns to the grain DataFrame.

| Column | Description |
|--------|-------------|
| `area_px` | Grain area in pixels |
| `area_um2` | Physical area (µm²) from step size |
| `centroid_x`, `centroid_y` | Centre of mass |
| `ellipse_a`, `ellipse_b` | Major/minor semi-axes |
| `ellipse_angle` | Orientation of major axis (degrees from +X) |
| `aspect_ratio` | `ellipse_a / ellipse_b` |
| `feret_max`, `feret_min` | Max/min caliper diameters |
| `convex_hull_area` | Area of convex hull |
| `solidity` | `area / convex_hull_area` |

## `santex.grains.symmetry` — Crystal symmetry

```python
from santex.grains.symmetry import get_laue_group
lg = get_laue_group("mmm")   # returns orix LaueGroup object
```

Supports all 11 Laue groups: `-1`, `2/m`, `mmm`, `4/m`, `4/mmm`, `-3`, `-3m`, `6/m`, `6/mmm`, `m-3`, `m-3m`.

## KAM — Kernel Average Misorientation

$$\text{KAM}(i) = \frac{1}{|\mathcal{N}(i)|} \sum_{j \in \mathcal{N}(i)} \omega(g_i, g_j)$$

where $\omega(g_i, g_j)$ is the misorientation angle between pixel $i$ and its neighbour $j$, and $\mathcal{N}(i)$ is the set of nearest-neighbour pixels in the same grain.

## GROD — Grain Reference Orientation Deviation

$$\text{GROD}(i) = \omega(g_i,\,\bar{g}_{\text{grain}})$$

where $\bar{g}_{\text{grain}}$ is the mean orientation of the grain containing pixel $i$.
