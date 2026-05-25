# `santex.ebsd` — EBSD Processing

```python
from santex.ebsd import EBSD
ebsd = EBSD("scan.ctf")
```

## Class `EBSD`

### Loading data

#### `__init__(filepath: str)`

Load an EBSD file. Supported formats: `.ctf`, `.ang`, `.osc`, `.txt`.

#### `get_ebsd_data() → pd.DataFrame`

Return the full data table. Columns include `Phase`, `X`, `Y`, `Euler1`, `Euler2`, `Euler3`, `MAD`, `BC`, `BS`.

```python
df = ebsd.get_ebsd_data()
print(df.head())
```

#### `get_ebsd_data_header() → dict`

Return header metadata: step size, grid dimensions, acquisition Euler angles, phase list.

#### `phases() → pd.DataFrame`

Return phase summary: index, name, fraction.

```python
phases = ebsd.phases()
# phase_index | phase_name | count | fraction
```

#### `phases_names() → pd.DataFrame`

Return just the phase-name column.

### Filtering

#### `filter_MAD(df: pd.DataFrame, threshold: float) → pd.DataFrame`

Remove pixels with MAD > threshold.

```python
df_clean = ebsd.filter_MAD(df, 0.8)
```

#### `filter_by_phase_number_list(df: pd.DataFrame, phase_list: list[int]) → pd.DataFrame`

Remove specific phase indices (e.g. `[0]` removes unindexed pixels).

#### `filter_by_grain_size(df: pd.DataFrame, phase_names: list[str], min_grain_size: int) → pd.DataFrame`

Remove grains smaller than `min_grain_size` pixels.

### Orientation operations

#### `get_euler_angles(phase: int, data: pd.DataFrame) → np.ndarray`

Return Euler angles (N × 3, degrees) for one phase.

```python
fo_euler = ebsd.get_euler_angles(phase=1, data=df_clean)
```

#### `rotate_ebsd(df: pd.DataFrame, angles: tuple) → pd.DataFrame`

Rotate the EBSD coordinate frame by Bunge ZXZ angles.

#### `plot_rotate_ebsd(sample_ref: list[str], ebsd_df: pd.DataFrame) → pd.DataFrame`

Rotate to a standard convention, e.g. `["x2east", "zOutOfPlane"]`.

### Grain reconstruction

#### `calc_grains(df, threshold, phase_names, downsampling_factor=1) → pd.DataFrame`

Segment grains by misorientation threshold. Returns the dataframe with a `grain_id` column added.

### Anisotropy from EBSD

#### `get_anisotropy_for_ebsd(cij_list, euler_list, density_list) → tuple[np.ndarray, float]`

Compute VRH-averaged elastic tensor and density for a multi-phase scan.

```python
average_cij, average_rho = ebsd.get_anisotropy_for_ebsd(
    cij_list=[cij_fo, cij_ens, cij_di],
    euler_list=[fo_euler, ens_euler, di_euler],
    density_list=[rho_fo, rho_ens, rho_di],
)
```

### Visualisation

#### `plot(data=None, cmap="viridis", save_image=False, image_filename=None, **kwargs)`

Plot the EBSD phase / Euler map.

#### `pf(df, phase, crystal_symmetry, uvw, ...)`

Pole figure for one phase and one crystallographic direction.

#### `ipf(df, phase, vector_sample, crystal_symmetry, ...)`

Inverse pole figure coloured by sample direction.

#### `ipf_colorkey(df, phase, crystal_symmetry)`

Show the IPF colour triangle key.
