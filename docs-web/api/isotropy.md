# `santex.isotropy` — Isotropic Velocities

Compute isotropic seismic velocities (*Vp*, *Vs*, *V*bulk) and elastic moduli at elevated *P*–*T* using the Hacker & Abers (2004) database.

```python
from santex.isotropy import Isotropy
iso = Isotropy()
```

## Class `Isotropy`

### `get_available_phases() → None`

Print all available minerals (60+ phases).

```python
iso.get_available_phases()
# Alpha-Quartz, Beta-Quartz, Forsterite, Diopside, ...
```

### `get_phase_constants(name: str) → dict`

Return the thermodynamic constants for one mineral.

```python
c = iso.get_phase_constants("Forsterite")
# {
#   'id': 'fo', 'name': 'Forsterite',
#   'rho0': 3222.0,        # kg/m³
#   'ao': 6.13e-5,         # thermal expansion coefficient
#   'akt0': 1.273e11,      # isothermal bulk modulus (Pa)
#   'dkdp': 4.2,           # dK/dP
#   'amu0': 8.16e10,       # shear modulus (Pa)
#   'dmudp': 1.6,          # dG/dP
#   'gam': 5.19,           # first Grüneisen parameter
#   'grun': 1.29,          # second Grüneisen parameter
#   'delt': 5.5,           # Debye temperature scaling
# }
```

### `calculate_seismic_properties(name, temperature, pressure, return_vp_vs_vbulk=True, return_aktout=True) → tuple`

Returns `(density, K_s, G, Vp, Vs, Vbulk, K_T)` at the given *T* (K) and *P* (GPa).

```python
rho, ks, g, vp, vs, vbulk, kt = iso.calculate_seismic_properties(
    "Forsterite", temperature=1500, pressure=2,
    return_vp_vs_vbulk=True, return_aktout=True,
)
print(f"Vp = {vp/1000:.3f} km/s  Vs = {vs/1000:.3f} km/s")
```

### Hashin–Shtrikman bounds

```python
phases = ['fo', 'di', 'en']
fractions = [0.62, 0.22, 0.16]

ph, fr = iso.set_modal_composition(phase_list=phases, fraction_list=fractions)

mid, upper, lower, rho = iso.hashin_shtrikman_bounds(
    phase_constant_list=ph,
    fraction_list=fr,
    temperature=1100,
    pressure=1.4,
    density_mix_calc=True,
)
# mid, upper, lower: dicts with keys 'K', 'G', 'vp', 'vs', 'vbulk', 'vp_vs_ratio'
```

## Physical model

Density at P, T:

$$\rho(P,T) = \rho_0 \exp\!\left[\int_0^P \frac{\mathrm{d}P'}{K_T(P')}\right]\exp\!\left[-\int_{300}^{T}\alpha(T')\,\mathrm{d}T'\right]$$

Seismic velocities:

$$V_p = \sqrt{\frac{K_s + \frac{4}{3}G}{\rho}}, \quad V_s = \sqrt{\frac{G}{\rho}}, \quad V_{\text{bulk}} = \sqrt{\frac{K_s}{\rho}}$$

Reference: Hacker & Abers (2004), *Geochem. Geophys. Geosyst.* 5.
