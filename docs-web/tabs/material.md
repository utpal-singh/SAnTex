# Material Database Tab

Browse the built-in mineral catalogue, set pressure and temperature conditions, and compute the single-crystal elastic tensor and anisotropy for any mineral.

## Panel layout

| Panel | Contents |
|-------|----------|
| **Left** | Mineral table (phase · crystal system · primary phase · density) |
| **Left** | P–T conditions spinboxes |
| **Left** | *Compute anisotropy* button + status label |
| **Right** | 6 × 6 Voigt matrix (GPa) |
| **Right** | Scalar anisotropy metrics table |
| **Right** | Plot Options + *Vp* stereonet |

## Controls

| Control | Range | Default | Description |
|---------|-------|---------|-------------|
| Pressure | 0 – 300 GPa | 0 | Confining pressure |
| Temperature | 0 – 3000 K | 300 | Absolute temperature |
| Mineral table | — | — | Click a row to select phase |

## Workflow

1. Click any mineral in the table (e.g. **Forsterite**).
2. Set **Pressure** and **Temperature** if different from ambient.
3. Click **Compute anisotropy**.
4. Read:
   - **Voigt matrix** — 6 × 6 elastic stiffness tensor in GPa.
   - **Anisotropy values** — scalar metrics (see table below).
   - ***Vp* stereonet** — P-wave velocity over the upper hemisphere.

## Anisotropy scalar metrics

| Metric | Symbol | Unit | Definition |
|--------|--------|------|------------|
| Max *Vp* | — | km/s | Maximum P-wave velocity |
| Min *Vp* | — | km/s | Minimum P-wave velocity |
| Max *Vs1* | — | km/s | Max fast shear velocity |
| Min *Vs1* | — | km/s | Min fast shear velocity |
| Max *Vs2* | — | km/s | Max slow shear velocity |
| Min *Vs2* | — | km/s | Min slow shear velocity |
| *P*-wave anisotropy | AVP | % | 200 × (Vp_max − Vp_min) / (Vp_max + Vp_min) |
| *S1*-wave anisotropy | AVS1 | % | 200 × (Vs1_max − Vs1_min) / (Vs1_max + Vs1_min) |
| *S2*-wave anisotropy | AVS2 | % | As above for Vs2 |
| AVS max | — | % | Max S-wave splitting |
| Vp / Vs1 | — | — | Velocity ratio |

## Pressure–temperature derivatives

The elastic constants are extrapolated linearly:

$$C_{ij}(p, T) = C_{ij}^0 + \frac{\partial C_{ij}}{\partial p}\,p + \frac{\partial C_{ij}}{\partial T}\,T$$

The database stores $C_{ij}^0$, $\partial C_{ij}/\partial p$, $\partial C_{ij}/\partial T$ for each mineral. Density is updated via:

$$\rho(p, T) = \rho_0 \left[1 + \frac{p}{K_0} - \alpha_0 (T - 300)\right]$$

::: warning Extrapolation
Derivatives are first-order. Large deviations from measured ranges (typically 0–10 GPa, 300–1500 K) reduce accuracy. Results outside calibrated ranges are labelled *extrapolated*.
:::

## Available minerals (excerpt)

| Phase | Crystal System | Primary group |
|-------|----------------|---------------|
| Forsterite | Orthorhombic | Olivine |
| Fayalite | Orthorhombic | Olivine |
| Enstatite | Orthorhombic | Olivine |
| Diopside | Monoclinic | Clinopyroxenes |
| Augite | Monoclinic | Clinopyroxenes |
| Hornblende (#1) | Monoclinic | Amphibole |
| Muscovite | Monoclinic | Phyllosilicates |
| Biotite | Monoclinic | Phyllosilicates |
| a_quartz_1 | Hexagonal/Trigonal | Quartz |
| Calcite_1 | Hexagonal/Trigonal | Calcite |
| Almandine-pyrope | Cubic | Garnet |
| Grossular | Cubic | Garnet |
| Albite | Triclinic | Plagioclase |
| An60 (Labradorite) | Triclinic | Plagioclase |

Full list: 55+ phases. Run `Material().available_phases()` in Python for the complete table.
