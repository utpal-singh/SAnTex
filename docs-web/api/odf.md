# `santex.odf` — ODF & Pole Figures

Orientation Distribution Function computation and texture analysis.

## Via `EBSD` class

```python
from santex.ebsd import EBSD
ebsd = EBSD("scan.ctf")
df = ebsd.get_ebsd_data()

# Pole figure
ebsd.pf(
    df=df,
    phase=1,
    crystal_symmetry='D2',     # or 'mmm', 'Ci', '-3m', etc.
    uvw=[0, 0, 1],             # Miller indices
    hemisphere='upper',
    axes_labels=['X', 'Y'],
)

# IPF
ebsd.ipf(
    df=df,
    phase=1,
    vector_sample=[0, 0, 1],   # sample reference direction
    crystal_symmetry='D2',
    vector_title='Z',
)
```

## `santex.odf.ODF` class (direct use)

```python
from santex.odf import ODF
import numpy as np

euler_angles = np.loadtxt("orientations.txt")   # N × 3, degrees

odf = ODF(euler_angles, crystal_symmetry='mmm', halfwidth_deg=5.0)
odf.compute()

print(f"Texture index J = {odf.texture_index():.2f}")
print(f"Entropy S       = {odf.entropy():.4f}")

# φ₂ section at φ₂ = 0°
odf.plot_phi2_section(phi2=0.0, colorscale='RdBu_r')

# Pole figure for [100]
odf.plot_pole_figure(hkl=[1, 0, 0])
```

## ODF theory

The de la Vallée Poussin KDE on SO(3):

$$f(g) = \frac{1}{N}\sum_{n=1}^{N} \kappa_\nu\!\left(\omega(g, g_n)\right)$$

where $\omega(g, g_n)$ is the disorientation angle and the kernel is:

$$\kappa_\nu(\omega) = \frac{\nu(\nu+1)}{4\pi(2\nu+1)} \left(\cos\tfrac{\omega}{2}\right)^{2\nu-2} \left(1 + 2\nu\cos\omega\right)$$

The halfwidth $h$ controls the kernel bandwidth: narrower halfwidth → sharper, noisier ODF; wider → smoother, less detail.

## Texture index

$$J = \int_{\text{SO}(3)} \!\! f(g)^2 \, \mathrm{d}g$$

| *J* value | Texture strength |
|-----------|-----------------|
| 1 | Perfect random (isotropic) |
| 2–5 | Weak |
| 5–15 | Moderate |
| >15 | Strong |
| ∞ | Single crystal |

## Export

```python
# POPLA epf (pole figure format)
odf.export_popla("texture.epf", phases=["Forsterite"])

# MTEX-compatible
odf.export_mtex("texture_mtex.txt")

# VPSC
odf.export_vpsc("vpsc_tex.txt", n_orientations=1000)
```
