# Anisotropy Tab

Compute seismic velocity surfaces and anisotropy metrics directly from a Voigt stiffness matrix. Visualise 2-D stereonets (Plotly) and a 3-D sphere (PyVista).

## Panel layout

| Panel | Contents |
|-------|----------|
| **Left** | 6 × 6 Voigt matrix input (editable cells) |
| **Left** | Density input, n θ / n φ resolution controls |
| **Left** | *Compute* button, scalar results label |
| **Right (top)** | 2-D stereonet (Plotly, switchable: Vp / Vs1 / Vs2 / AVS) |
| **Right (bottom)** | 3-D velocity surface (PyVista) |

## Controls

| Control | Default | Description |
|---------|---------|-------------|
| Voigt matrix cells | 0 | C₁₁…C₆₆ in GPa — paste or type values |
| Density | 0 kg/m³ | Used for velocity: V = √(C/ρ) |
| n θ | 90 | Number of polar angle samples (latitude) |
| n φ | 180 | Number of azimuthal angle samples (longitude) |
| Plot quantity | Vp | Dropdown: Vp, Vs1, Vs2, AVS, Vp/Vs1 |

## Velocity computation

Velocities are obtained by solving the **Christoffel equation** for each wave-propagation direction **n̂**:

$$\det\!\left[C_{ikjl}\,n_k n_l - \rho V^2\,\delta_{ij}\right] = 0$$

The three real eigenvalues of the Christoffel matrix $M_{ij} = C_{ikjl}\,n_k n_l$ give:

- $V_p^2 = \lambda_{\max} / \rho$ — fastest eigenvalue → P wave
- $V_{s1}^2 = \lambda_{\text{mid}} / \rho$ — middle eigenvalue → fast S wave
- $V_{s2}^2 = \lambda_{\min} / \rho$ — slowest eigenvalue → slow S wave

## Workflow

1. Fill the 6 × 6 Voigt matrix (you can copy from the Material Database tab after computing).
2. Enter the density (kg/m³).
3. Adjust **n θ** and **n φ** — higher values give smoother surfaces but slower computation.
4. Click **Compute**.
5. Switch between **Vp / Vs1 / Vs2 / AVS** in the dropdown to explore different surfaces.
6. The 3-D panel updates simultaneously.

## Resolution guide

| n θ × n φ | Points | Speed | Use case |
|-----------|--------|-------|----------|
| 45 × 90 | 4 050 | < 0.1 s | Preview |
| 90 × 180 | 16 200 | ~0.5 s | Standard |
| 180 × 360 | 64 800 | ~3 s | Publication |

::: tip Tip
Right-click the 3-D surface to change colour map, add lighting, or export a screenshot. Use the toolbar's SVG button on the 2-D stereonet for publication-quality exports.
:::
