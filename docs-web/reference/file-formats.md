# File Formats

## Input formats

| Format | Extension | Source software |
|--------|-----------|-----------------|
| HKL Channel 5 | `.ctf` | Oxford Instruments AZtec / HKL |
| TSL / EDAX | `.ang` | EDAX OIM Analysis |
| Oxford binary | `.osc` | Oxford Instruments |
| Generic text | `.txt` | Tab-separated with header |

### `.ctf` structure

```
Channel Text File
Prj     my_project
Author  J. Doe
JobMode Grid
XCells  400
YCells  300
XStep   0.1
YStep   0.1
AcqE1   0
AcqE2   0
AcqE3   0
Euler angles refer to Sample Coordinate system (CS0)!
Phases 2
3.956;3.956;3.956   90;90;90   Spinel ...
4.756;10.227;5.994  90;105.1;90  Monoclinic ...
Phase X    Y    Bands Error Euler1  Euler2  Euler3  MAD  BC  BS
1     0.1  0.1  9     0     125.3   58.2    312.1   0.41 188 255
1     0.2  0.1  8     0     124.9   57.8    311.7   0.38 190 252
...
```

Key columns:

| Column | Description |
|--------|-------------|
| `Phase` | Integer phase index (0 = unindexed) |
| `X`, `Y` | Pixel coordinates (µm) |
| `Euler1`, `Euler2`, `Euler3` | Bunge ZXZ Euler angles φ₁, Φ, φ₂ (degrees) |
| `MAD` | Mean Angular Deviation — indexing quality (lower = better) |
| `BC` | Band Contrast — diffraction pattern quality |
| `BS` | Band Slope |

### `.ang` structure

TSL format columns: `phi1 PHI phi2 x y IQ CI phase_id detector_intensity fit`

- `IQ` = Image Quality (band contrast equivalent)
- `CI` = Confidence Index (MAD equivalent, 0–1)

## Output formats

| Format | Extension | Tab | Consumer |
|--------|-----------|-----|---------|
| POPLA pole figures | `.epf` | ODF & PF | EPSC, popLA |
| MTEX text | `.txt` | ODF & PF | MTEX Matlab toolbox |
| VPSC orientations | `.txt` | ODF & PF | Los Alamos VPSC |
| Grain summary | `.csv` | Grains | Excel, Python, R |
| Neper ORI | `.ori` | Grains | Neper polycrystal |
| Neper TESS | `.tess` | Grains | Neper FEM |
| HKL CTF | `.ctf` | EBSD / Grains | AZtec, MTEX, HKL |

### POPLA `.epf` structure

Pole figure data format used by the EPSC and popLA codes:

```
pole figure for Forsterite   phase 1
100
  90.00   0.00   1.0000
  90.00  10.00   0.9832
  ...
```

Each line: θ (degrees), φ (degrees), intensity.

### VPSC orientations

```
  500  weighted   # number of orientations
  125.3   58.2  312.1   1.0   # phi1 PHI phi2 weight
  124.9   57.8  311.7   1.0
  ...
```
