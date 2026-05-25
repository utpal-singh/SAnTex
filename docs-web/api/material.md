# `santex.material` — Material Database

```python
from santex.material import Material
mat = Material()
```

## Class `Material`

### `available_phases() → pd.DataFrame`

Return the full mineral catalogue.

```python
df = mat.available_phases()
print(df.to_string())
```

Columns: `phase`, `crystal_system`, `primary_phase`, `density`.

### `get_properties_by_phase(phase: str) → dict`

Return all stored properties for one mineral.

```python
props = mat.get_properties_by_phase("Diopside")
# Keys: Crystal System, Primary Phase, Phase, Density(g/cm3),
#       C11..C66, Crystal Reference Frame, Study
```

### `get_voigt_matrix(phase: str) → np.ndarray`

Return the 6 × 6 Voigt stiffness matrix (GPa) at ambient conditions.

```python
cij = mat.get_voigt_matrix("Forsterite")   # shape (6, 6)
```

### `voigt_high_PT(phase: str, PRESSURE: float, TEMP: float) → np.ndarray`

Return the 6 × 6 matrix extrapolated to pressure (GPa) and temperature (°C or K — check database notes).

```python
cij = mat.voigt_high_PT("Forsterite", PRESSURE=2, TEMP=1500)
```

### `load_density(phase: str, pressure: float = 0, temperature: float = 300) → float`

Return density (kg/m³) at the given P–T.

```python
rho = mat.load_density("Forsterite", pressure=2, temperature=1500)
# → 3301.4  (kg/m³)
```

### `modal_rock(phases: list[str], fractions: list[float], pressure: float, temperature: float) → tuple[np.ndarray, float]`

Compute the Voigt-average elastic matrix and density for a mixture.

```python
cij, rho = mat.modal_rock(
    ["Forsterite", "Diopside", "Enstatite"],
    [0.65, 0.20, 0.15],
    pressure=2,
    temperature=1000,
)
# cij shape (6, 6) in GPa units (multiply by 1e9 for Pa)
```

### `get_density(phase: str) → float | None`

Ambient density in g/cm³ (multiply by 1000 for kg/m³).
