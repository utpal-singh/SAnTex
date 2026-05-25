# Modal Rock Tab

Compute the elastic tensor of a rock from its modal mineral composition (volume fractions) without EBSD data. Assumes a **random texture** (isotropic orientation distribution).

## Panel layout

| Panel | Contents |
|-------|----------|
| **Left top** | Phase/fraction input table |
| **Left** | Add / Remove row buttons |
| **Left** | P, T, averaging method |
| **Left** | *Compute Modal Rock* button |
| **Right top** | Voigt matrix + scalar metrics |
| **Right** | Hashin–Shtrikman sub-tab |
| **Right** | Plotly *Vp* stereonet |

## Controls

| Control | Default | Description |
|---------|---------|-------------|
| Phase name | — | Must match a name in the Material database |
| Volume fraction | — | 0 – 100; all rows must sum to 100 |
| Pressure | 0 GPa | |
| Temperature | 300 K | |
| Averaging | Hill | Voigt / Reuss / Hill |

## Workflow

1. Click **Add phase** and type a mineral name (autocomplete from database).
2. Enter its volume fraction (%).
3. Repeat until all phases are added and fractions sum to 100 %.
4. Set P and T.
5. Click **Compute Modal Rock**.
6. The Voigt matrix, anisotropy values, and *Vp* stereonet appear.
7. Switch to the **Hashin–Shtrikman** sub-tab for isotropic bounds.

::: tip Example
A typical upper mantle peridotite:
| Phase | Fraction |
|-------|----------|
| Forsterite | 65 % |
| Enstatite | 20 % |
| Diopside | 10 % |
| Pyrope | 5 % |
:::

## Hashin–Shtrikman bounds

The H–S bounds give the tightest possible isotropic stiffness range for a random mixture of *n* phases without assumptions about geometry:

$$K^{\text{HS}-} \leq K_{\text{rock}} \leq K^{\text{HS}+}$$
$$G^{\text{HS}-} \leq G_{\text{rock}} \leq G^{\text{HS}+}$$

From *K* and *G*, isotropic velocities follow:

$$V_p = \sqrt{\frac{K + \tfrac{4}{3}G}{\rho}}, \qquad V_s = \sqrt{\frac{G}{\rho}}$$

The sub-tab shows:
- Upper H–S bound ($K^+$, $G^+$, $V_p^+$, $V_s^+$)
- Lower H–S bound ($K^-$, $G^-$, $V_p^-$, $V_s^-$)
- VRH mean
- Density ρ (weighted average)

::: info
H–S bounds use the Hacker & Abers (2004) isotropic velocity database, which is separate from the anisotropic Voigt matrix database used for the main computation.
:::
