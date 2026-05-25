# `santex.anisotropy` — Seismic Anisotropy

```python
from santex.anisotropy import Anisotropy
anis = Anisotropy(cij_pa, rho)   # Pa, kg/m³
```

## Class `Anisotropy`

### `__init__(cij: np.ndarray, rho: float)`

Initialise with:
- `cij` — 6 × 6 Voigt stiffness matrix **in Pa** (not GPa)
- `rho` — density in kg/m³

```python
from santex.material import Material
from santex.anisotropy import Anisotropy

mat = Material()
cij = mat.voigt_high_PT("Forsterite", PRESSURE=2, TEMP=1500)
rho = mat.load_density("Forsterite", 2, 1500)

anis = Anisotropy(cij * 1e9, rho)
```

### `anisotropy_values() → dict`

Compute and print all scalar anisotropy metrics. Returns a dict with keys:

| Key | Unit | Description |
|-----|------|-------------|
| `maxvp` | m/s | Maximum P-wave velocity |
| `minvp` | m/s | Minimum P-wave velocity |
| `maxvs1` | m/s | Maximum fast S-wave velocity |
| `minvs1` | m/s | Minimum fast S-wave velocity |
| `maxvs2` | m/s | Maximum slow S-wave velocity |
| `minvs2` | m/s | Minimum slow S-wave velocity |
| `max_vs_anisotropy_percent` | % | Maximum shear-wave splitting |
| `min_vs_anisotropy_percent` | % | Minimum shear-wave splitting |
| `p_wave_anisotropy_percent` | % | AVP = 200(Vp_max−Vp_min)/(Vp_max+Vp_min) |
| `s1_wave_anisotropy_percent` | % | AVS1 |
| `s2_wave_anisotropy_percent` | % | AVS2 |
| `maxdvs` | m/s | Max S-wave velocity difference |
| `AVpVs1` | — | Max Vp/Vs1 ratio |

```python
vals = anis.anisotropy_values()
print(f"AVP = {vals['p_wave_anisotropy_percent']:.2f} %")
```

### `plot(colormap="RdBu_r", step=180, savefig=False, figname=None, dpi=300, save_format="svg")`

Matplotlib 2-D stereonet (four panels: Vp, Vs1, Vs2, AVS).

### `plotter_vp(rho, cij)` / `plotter_vs1(...)` / `plotter_vs2(...)` / `plotter_vs_splitting(...)`

PyVista 3-D interactive surface for each wave type.

```python
anis.plotter_vp(rho, cij)     # opens interactive 3-D window
```

## Velocity formula

Each velocity is derived from the eigenvalues λ of the Christoffel matrix:

$$M_{ij}(\hat{n}) = C_{ikjl}\,n_k\,n_l$$

$$\det(M_{ij} - \rho V^2 \delta_{ij}) = 0 \implies V = \sqrt{\lambda / \rho}$$
