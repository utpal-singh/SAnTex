"""
Pure-computation EBSD operations (no UI, no Plotly/matplotlib).

All functions take / return NumPy arrays or pandas DataFrames so they
can be called from any frontend.

Functions
---------
orix_symmetry(name)            → orix Symmetry object
ipf_map_colors(data, ...)      → dict(x, y, rgb_hex, r, g, b)
bc_bs_map_data(data)           → dict(x, y, bc, bs)
compute_kam(data, ...)         → (N,) array of KAM angles (degrees)
compute_mis2mean(data, phase)  → (N,) array of Mis2Mean angles (degrees)
denoise_mean_filter(data, ...) → updated DataFrame
fill_missing_data(data, ...)   → updated DataFrame
regrid(data, target_step)      → new regular-grid DataFrame
select_by_condition(data, ...) → filtered DataFrame
get_line_profile(data, ...)    → dict with profile scalars
simulate_ebsd(...)             → synthetic EBSD DataFrame
"""

from __future__ import annotations
import numpy as np
import pandas as pd


# ---------------------------------------------------------------------------
# Symmetry helper — avoids the 300-line if/elif in crystal.py
# ---------------------------------------------------------------------------

_ORIX_SYM_ALIASES: dict[str, str] = {
    # triclinic
    "1": "C1", "Ci": "Ci", "-1": "Ci",
    # monoclinic
    "2": "C2", "m": "Cs", "2/m": "C2h",
    # orthorhombic
    "222": "D2", "mm2": "C2v", "mmm": "D2h",
    # tetragonal
    "4": "C4", "-4": "S4", "4/m": "C4h",
    "422": "D4", "4mm": "C4v", "-42m": "D2d", "-4m2": "D2d", "4/mmm": "D4h",
    # trigonal
    "3": "C3", "-3": "C3i",
    "32": "D3", "3m": "C3v", "-3m": "D3d",
    "321": "D3", "3m1": "C3v", "-3m1": "D3d",
    "312": "D3", "31m": "C3v", "-31m": "D3d",
    # hexagonal
    "6": "C6", "-6": "C3h", "6/m": "C6h",
    "622": "D6", "6mm": "C6v", "-62m": "D3h", "-6m2": "D3h", "6/mmm": "D6h",
    # cubic
    "23": "T", "m-3": "Th",
    "432": "O", "-43m": "Td", "m-3m": "Oh",
}


def orix_symmetry(name: str):
    """Return the orix symmetry object for *name*.

    Accepts both the canonical class names (``'D2h'``, ``'Oh'``, …) and the
    Hermann–Mauguin or Schoenflies short forms (``'mmm'``, ``'m-3m'``, …).
    Falls back to ``D2h`` (mmm) if the name is unrecognised.
    """
    from orix.quaternion import symmetry as sym

    _MAP = {
        "C1": sym.C1, "Ci": sym.Ci,
        "C2": sym.C2, "Cs": sym.Cs, "C2h": sym.C2h,
        "D2": sym.D2, "C2v": sym.C2v, "D2h": sym.D2h,
        "C3": sym.C3, "C3i": sym.C3i,
        "D3": sym.D3, "C3v": sym.C3v, "D3d": sym.D3d,
        "C4": sym.C4, "S4": sym.S4, "C4h": sym.C4h,
        "D4": sym.D4, "C4v": sym.C4v, "D2d": sym.D2d, "D4h": sym.D4h,
        "C6": sym.C6, "C3h": sym.C3h, "C6h": sym.C6h,
        "D6": sym.D6, "C6v": sym.C6v, "D3h": sym.D3h, "D6h": sym.D6h,
        "T": sym.T, "Th": sym.Th,
        "O": sym.O, "Td": sym.Td, "Oh": sym.Oh,
    }
    key = _ORIX_SYM_ALIASES.get(name, name)
    return _MAP.get(key, sym.D2h)


# ---------------------------------------------------------------------------
# IPF map colours
# ---------------------------------------------------------------------------

def ipf_map_colors(
    data: pd.DataFrame,
    phase_index: int = 1,
    direction: str = "ND",
    crystal_symmetry: str = "D2h",
) -> dict | None:
    """
    Compute per-pixel IPF colour for a single phase using orix.

    Parameters
    ----------
    data            : EBSD DataFrame (Phase, X, Y, Euler1, Euler2, Euler3)
    phase_index     : phase number (>0)
    direction       : ``'ND'`` (Z), ``'RD'`` (X), or ``'TD'`` (Y)
    crystal_symmetry: e.g. ``'D2h'``, ``'Oh'``, ``'mmm'``, ``'m-3m'``

    Returns
    -------
    dict with keys x, y, r, g, b (float 0-1), rgb_hex (list of str),
    or ``None`` if the phase has no indexed pixels.
    """
    from orix.quaternion import Orientation
    from orix.vector import Vector3d
    from orix import plot as orix_plot

    _DIR = {"ND": Vector3d([0, 0, 1]), "RD": Vector3d([1, 0, 0]),
            "TD": Vector3d([0, 1, 0])}
    vec = _DIR.get(direction.upper(), _DIR["ND"])

    phase_data = data[data["Phase"] == phase_index].copy()
    if len(phase_data) == 0:
        return None

    euler = phase_data[["Euler1", "Euler2", "Euler3"]].to_numpy()
    sym_obj = orix_symmetry(crystal_symmetry)
    ori = Orientation.from_euler(np.deg2rad(euler), sym_obj)

    ipfkey = orix_plot.IPFColorKeyTSL(sym_obj, direction=vec)
    rgb = ipfkey.orientation2color(ori)          # (N, 3) float in [0, 1]
    rgb = np.clip(rgb, 0.0, 1.0)

    hex_colors = [
        "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))
        for r, g, b in rgb
    ]
    return {
        "x":       phase_data["X"].values,
        "y":       phase_data["Y"].values,
        "r":       rgb[:, 0],
        "g":       rgb[:, 1],
        "b":       rgb[:, 2],
        "rgb_hex": hex_colors,
    }


# ---------------------------------------------------------------------------
# BC / BS map
# ---------------------------------------------------------------------------

def bc_bs_map_data(data: pd.DataFrame) -> dict | None:
    """Return X, Y, BC, BS arrays for all pixels (or those with Phase > 0)."""
    if data is None or "BC" not in data.columns:
        return None
    mask = data["Phase"] > 0
    d = data[mask]
    result = {"x": d["X"].values, "y": d["Y"].values}
    if "BC" in d.columns:
        result["bc"] = d["BC"].values.astype(float)
    if "BS" in d.columns:
        result["bs"] = d["BS"].values.astype(float)
    return result


# ---------------------------------------------------------------------------
# KAM — from raw EBSD data (no grain reconstruction required)
# ---------------------------------------------------------------------------

def compute_kam(
    data: pd.DataFrame,
    kernel_size: int = 1,
    max_angle_deg: float = 5.0,
    same_phase_only: bool = True,
) -> np.ndarray:
    """
    Kernel Average Misorientation per pixel, computed directly from EBSD data.

    Builds a cKDTree on (X, Y) positions and searches within
    ``kernel_size × step`` radius for each pixel.

    Parameters
    ----------
    data            : EBSD DataFrame
    kernel_size     : neighbourhood radius in units of the grid step
    max_angle_deg   : misorientation cut-off; pairs above this are ignored
    same_phase_only : only compare pixels of the same phase

    Returns
    -------
    (N,) float array of KAM angles in degrees; NaN for phase-0 pixels.
    """
    from scipy.spatial import cKDTree

    n = len(data)
    x = data["X"].values.astype(float)
    y = data["Y"].values.astype(float)
    phase = data["Phase"].values.astype(int)

    phi1 = np.radians(data["Euler1"].values)
    PHI  = np.radians(data["Euler2"].values)
    phi2 = np.radians(data["Euler3"].values)
    R    = _euler_to_rotation_batch(phi1, PHI, phi2)   # (N,3,3)

    # Estimate step size from median of nearest neighbour distances
    step = _estimate_step(x, y)
    radius = (kernel_size + 0.5) * step

    tree = cKDTree(np.column_stack([x, y]))
    kam     = np.full(n, np.nan)
    kam_sum = np.zeros(n)
    kam_cnt = np.zeros(n, dtype=np.int32)

    for i, (xi, yi) in enumerate(zip(x, y)):
        if phase[i] == 0:
            continue
        idx = tree.query_ball_point([xi, yi], radius)
        idx = [j for j in idx if j != i and phase[j] > 0]
        if same_phase_only:
            idx = [j for j in idx if phase[j] == phase[i]]
        if not idx:
            kam[i] = 0.0
            continue
        M   = np.einsum('ij,kjl->kil', R[i].T, R[np.array(idx)])
        tr  = np.clip(np.einsum('kii->k', M), -1.0, 3.0)
        ang = np.degrees(np.arccos((tr - 1.0) / 2.0))
        valid = ang < max_angle_deg
        if valid.any():
            kam[i] = ang[valid].mean()
        else:
            kam[i] = 0.0

    return kam


# ---------------------------------------------------------------------------
# Mis2Mean — misorientation from phase / region mean
# ---------------------------------------------------------------------------

def compute_mis2mean(
    data: pd.DataFrame,
    phase_index: int = 1,
) -> np.ndarray:
    """
    Per-pixel misorientation from the mean orientation of *phase_index*.

    The mean orientation is computed as the dominant eigenvector of the
    sum of outer products of unit quaternions (Markley et al. method).

    Returns
    -------
    (N,) float array; values for phase_index pixels are Mis2Mean angles in
    degrees; all other pixels get NaN.
    """
    n = len(data)
    phase = data["Phase"].values.astype(int)
    mask = phase == phase_index
    result = np.full(n, np.nan)

    if not mask.any():
        return result

    phi1 = np.radians(data["Euler1"].values[mask])
    PHI  = np.radians(data["Euler2"].values[mask])
    phi2 = np.radians(data["Euler3"].values[mask])
    R = _euler_to_rotation_batch(phi1, PHI, phi2)

    # Quaternion averaging via eigenvalue method
    Q = _rotation_to_quaternion_batch(R)           # (K,4)
    M = Q.T @ Q
    _, vecs = np.linalg.eigh(M)
    q_mean = vecs[:, -1]                           # dominant eigenvector = mean
    R_mean = _quaternion_to_rotation(q_mean)       # (3,3)

    # Misorientation from mean
    Mdiff = np.einsum('ij,kjl->kil', R_mean.T, R)
    tr    = np.clip(np.einsum('kii->k', Mdiff), -1.0, 3.0)
    angles = np.degrees(np.arccos((tr - 1.0) / 2.0))
    result[mask] = angles
    return result


# ---------------------------------------------------------------------------
# Denoising — moving-average filter in orientation space
# ---------------------------------------------------------------------------

def denoise_mean_filter(
    data: pd.DataFrame,
    kernel_size: int = 1,
    max_angle_deg: float = 5.0,
    fill_zero_phase: bool = False,
) -> pd.DataFrame:
    """
    Replace each pixel's orientation with the quaternion mean of its
    same-phase neighbours within ``kernel_size × step`` radius, provided
    the angular spread is < ``max_angle_deg``.

    Phase-0 pixels are left as-is unless *fill_zero_phase* is True.

    Returns a **copy** of *data* with updated Euler1/2/3 columns.
    """
    from scipy.spatial import cKDTree

    d = data.copy()
    x = d["X"].values.astype(float)
    y = d["Y"].values.astype(float)
    phase = d["Phase"].values.astype(int)
    step = _estimate_step(x, y)
    radius = (kernel_size + 0.5) * step

    phi1 = np.radians(d["Euler1"].values)
    PHI  = np.radians(d["Euler2"].values)
    phi2 = np.radians(d["Euler3"].values)
    R = _euler_to_rotation_batch(phi1, PHI, phi2)
    Q = _rotation_to_quaternion_batch(R)           # (N,4)

    tree = cKDTree(np.column_stack([x, y]))
    new_e1 = d["Euler1"].values.copy()
    new_e  = d["Euler2"].values.copy()
    new_e2 = d["Euler3"].values.copy()

    for i in range(len(d)):
        if phase[i] == 0 and not fill_zero_phase:
            continue
        idx = tree.query_ball_point([x[i], y[i]], radius)
        idx = [j for j in idx if phase[j] == phase[i] or
               (fill_zero_phase and phase[j] > 0)]
        if len(idx) < 2:
            continue
        # Quaternion mean
        Qn = Q[np.array(idx)]
        # Ensure consistent hemisphere
        if Q[i] @ Qn[0] < 0:
            Qn = -Qn
        M_outer = Qn.T @ Qn
        _, vecs = np.linalg.eigh(M_outer)
        q_mean = vecs[:, -1]
        # Check angular spread before accepting
        angles = np.degrees(np.arccos(np.clip(np.abs(Qn @ q_mean), 0, 1)) * 2)
        if angles.mean() > max_angle_deg:
            continue
        e1, e, e2 = _rotation_to_euler(_quaternion_to_rotation(q_mean))
        new_e1[i] = np.degrees(e1)
        new_e[i]  = np.degrees(e)
        new_e2[i] = np.degrees(e2)

    d["Euler1"] = new_e1
    d["Euler2"] = new_e
    d["Euler3"] = new_e2
    return d


# ---------------------------------------------------------------------------
# Fill missing data
# ---------------------------------------------------------------------------

def fill_missing_data(
    data: pd.DataFrame,
    method: str = "nearest",
) -> pd.DataFrame:
    """
    Fill phase-0 (unindexed) pixels by interpolation from indexed neighbours.

    Parameters
    ----------
    method : ``'nearest'`` — copy orientation of the spatially nearest
             indexed pixel (fast);
             ``'mean'`` — quaternion mean of the *n* nearest indexed pixels
             (smoother but slower).

    Returns a **copy** of *data*.
    """
    from scipy.spatial import cKDTree

    d = data.copy()
    phase = d["Phase"].values.astype(int)
    x = d["X"].values.astype(float)
    y = d["Y"].values.astype(float)

    indexed_mask = phase > 0
    unindexed_mask = phase == 0

    if not unindexed_mask.any():
        return d   # nothing to fill

    indexed_idx = np.where(indexed_mask)[0]
    tree = cKDTree(np.column_stack([x[indexed_idx], y[indexed_idx]]))

    q_x = x[unindexed_mask]
    q_y = y[unindexed_mask]

    if method == "nearest":
        _, nn = tree.query(np.column_stack([q_x, q_y]), k=1)
        src = indexed_idx[nn]
        for col in ["Euler1", "Euler2", "Euler3", "Phase"]:
            if col in d.columns:
                d.loc[unindexed_mask, col] = d[col].values[src]
    elif method == "mean":
        k = min(5, len(indexed_idx))
        _, nns = tree.query(np.column_stack([q_x, q_y]), k=k)
        src = indexed_idx[nns]   # (M, k)

        phi1_s = np.radians(d["Euler1"].values[src])
        PHI_s  = np.radians(d["Euler2"].values[src])
        phi2_s = np.radians(d["Euler3"].values[src])
        R_s = _euler_to_rotation_batch(
            phi1_s.ravel(), PHI_s.ravel(), phi2_s.ravel()
        ).reshape(len(q_x), k, 3, 3)
        Q_s = _rotation_to_quaternion_batch(
            R_s.reshape(-1, 3, 3)
        ).reshape(len(q_x), k, 4)

        un_idx = np.where(unindexed_mask)[0]
        for i in range(len(un_idx)):
            Qn = Q_s[i]
            M_outer = Qn.T @ Qn
            _, vecs = np.linalg.eigh(M_outer)
            q_mean = vecs[:, -1]
            e1, e, e2 = _rotation_to_euler(_quaternion_to_rotation(q_mean))
            d.at[un_idx[i], "Euler1"] = np.degrees(e1)
            d.at[un_idx[i], "Euler2"] = np.degrees(e)
            d.at[un_idx[i], "Euler3"] = np.degrees(e2)
            # Assign most-common phase among neighbours
            ph_vals = d["Phase"].values[src[i]]
            d.at[un_idx[i], "Phase"] = int(np.bincount(ph_vals[ph_vals > 0]).argmax()
                                           if (ph_vals > 0).any() else 1)
    return d


# ---------------------------------------------------------------------------
# Regridding / interpolation
# ---------------------------------------------------------------------------

def regrid(
    data: pd.DataFrame,
    target_step: float | None = None,
    scale: float | None = None,
) -> pd.DataFrame:
    """
    Resample EBSD data onto a new regular Cartesian grid.

    Either *target_step* (µm, absolute) or *scale* (multiplier) must be given.
    E.g. ``scale=0.5`` gives twice the resolution (interpolated); ``scale=2``
    coarsens by ×2.

    For each new grid point the orientation is set to the nearest original
    pixel (no interpolation in orientation space — nearest-neighbour).
    Scalar properties (BC, BS, MAD) are bilinearly interpolated on the
    original grid.

    Returns a **new** DataFrame on the regular grid.
    """
    from scipy.spatial import cKDTree
    from scipy.interpolate import griddata

    x = data["X"].values.astype(float)
    y = data["Y"].values.astype(float)

    step = _estimate_step(x, y)
    if target_step is not None:
        new_step = float(target_step)
    elif scale is not None:
        new_step = step * float(scale)
    else:
        return data.copy()

    xi = np.arange(x.min(), x.max() + new_step / 2, new_step)
    yi = np.arange(y.min(), y.max() + new_step / 2, new_step)
    XI, YI = np.meshgrid(xi, yi)
    xy_new = np.column_stack([XI.ravel(), YI.ravel()])

    tree = cKDTree(np.column_stack([x, y]))
    dists, nn = tree.query(xy_new, k=1)

    # Discard new points too far from any real pixel (>1.5× new_step)
    valid = dists < 1.5 * new_step
    nn_v  = nn[valid]

    rows = data.iloc[nn_v].copy().reset_index(drop=True)
    rows["X"] = XI.ravel()[valid]
    rows["Y"] = YI.ravel()[valid]

    # Optionally smooth scalar columns via griddata
    for col in ["BC", "BS", "MAD"]:
        if col in data.columns:
            try:
                vals = griddata(
                    np.column_stack([x, y]),
                    data[col].values.astype(float),
                    xy_new[valid],
                    method="nearest",
                )
                rows[col] = vals
            except Exception:
                pass

    return rows


# ---------------------------------------------------------------------------
# Select by condition / index
# ---------------------------------------------------------------------------

def select_by_condition(
    data: pd.DataFrame,
    phase: int | list[int] | None = None,
    mad_max: float | None = None,
    bc_min: float | None = None,
    bs_min: float | None = None,
    index_list: list[int] | np.ndarray | None = None,
    x_range: tuple[float, float] | None = None,
    y_range: tuple[float, float] | None = None,
) -> pd.DataFrame:
    """
    Return a filtered **view** of *data*.

    All conditions are combined with AND logic.
    """
    mask = pd.Series(True, index=data.index)

    if phase is not None:
        phases = [phase] if isinstance(phase, int) else list(phase)
        mask &= data["Phase"].isin(phases)
    if mad_max is not None and "MAD" in data.columns:
        mask &= data["MAD"] <= mad_max
    if bc_min is not None and "BC" in data.columns:
        mask &= data["BC"] >= bc_min
    if bs_min is not None and "BS" in data.columns:
        mask &= data["BS"] >= bs_min
    if index_list is not None:
        idx_set = set(index_list)
        mask &= data.index.isin(idx_set)
    if x_range is not None:
        mask &= (data["X"] >= x_range[0]) & (data["X"] <= x_range[1])
    if y_range is not None:
        mask &= (data["Y"] >= y_range[0]) & (data["Y"] <= y_range[1])

    return data[mask].copy().reset_index(drop=True)


# ---------------------------------------------------------------------------
# Line profiles
# ---------------------------------------------------------------------------

def get_line_profile(
    data: pd.DataFrame,
    x0: float, y0: float,
    x1: float, y1: float,
    n_points: int = 100,
    scalars: list[str] | None = None,
) -> dict:
    """
    Extract a line profile through the EBSD map from (x0,y0) to (x1,y1).

    For each of *n_points* equidistant sample positions along the line the
    nearest pixel value is returned.

    Returns
    -------
    dict with keys:
        ``distance``  : (n_points,) cumulative distance from start (µm)
        ``x``, ``y``  : (n_points,) sample coordinates
        *plus one array per scalar column requested*
    """
    from scipy.spatial import cKDTree

    t     = np.linspace(0, 1, n_points)
    xs    = x0 + t * (x1 - x0)
    ys    = y0 + t * (y1 - y0)
    dist  = t * np.sqrt((x1 - x0)**2 + (y1 - y0)**2)

    tree = cKDTree(np.column_stack([data["X"].values, data["Y"].values]))
    _, nn = tree.query(np.column_stack([xs, ys]), k=1)

    if scalars is None:
        scalars = [c for c in data.columns if c not in ("X", "Y")]

    result: dict = {"distance": dist, "x": xs, "y": ys}
    for col in scalars:
        if col in data.columns:
            result[col] = data[col].values[nn]
    return result


# ---------------------------------------------------------------------------
# Orientation scatter (Euler space 2D / 3D)
# ---------------------------------------------------------------------------

def orientation_scatter_data(
    data: pd.DataFrame,
    phase_index: int = 1,
    axes: tuple[str, str] = ("Euler1", "Euler2"),
    max_points: int = 5000,
) -> dict | None:
    """Return up to *max_points* orientation scatter data for the chosen axes."""
    d = data[data["Phase"] == phase_index]
    if len(d) == 0:
        return None
    if len(d) > max_points:
        d = d.sample(max_points, random_state=0)
    ax1, ax2 = axes
    result: dict = {}
    for ax in (ax1, ax2, "Euler1", "Euler2", "Euler3"):
        if ax in d.columns:
            result[ax] = d[ax].values
    return result


# ---------------------------------------------------------------------------
# Simulation — generate synthetic EBSD DataFrame
# ---------------------------------------------------------------------------

def simulate_ebsd(
    euler_angles: np.ndarray,
    n_cols: int = 50,
    n_rows: int = 50,
    step: float = 1.0,
    phase_id: int = 1,
    noise_deg: float = 0.5,
    rng_seed: int = 42,
) -> pd.DataFrame:
    """
    Build a synthetic EBSD dataset from a set of prototype orientations.

    Each grid pixel is assigned one of the orientations in *euler_angles*
    (chosen at random), optionally perturbed by *noise_deg* degrees.

    Parameters
    ----------
    euler_angles : (M, 3) Bunge ZXZ angles in degrees
    n_cols, n_rows : grid size in pixels
    step          : pixel size in µm
    phase_id      : phase index to assign to all pixels
    noise_deg     : isotropic orientation noise (°)
    rng_seed      : random seed for reproducibility
    """
    rng = np.random.default_rng(rng_seed)
    n_total = n_cols * n_rows

    # Draw with replacement from prototype orientations
    idx = rng.integers(0, len(euler_angles), size=n_total)
    e = euler_angles[idx].astype(float)

    # Add Gaussian noise
    if noise_deg > 0:
        noise = rng.normal(0, noise_deg, size=(n_total, 3))
        e += noise
        e[:, 0] %= 360.0
        e[:, 1]  = np.clip(e[:, 1], 0.0, 90.0)
        e[:, 2] %= 360.0

    # Regular grid
    xi = np.arange(n_cols) * step
    yi = np.arange(n_rows) * step
    XI, YI = np.meshgrid(xi, yi)

    df = pd.DataFrame({
        "Phase":  phase_id,
        "X":      XI.ravel(),
        "Y":      YI.ravel(),
        "Euler1": e[:, 0],
        "Euler2": e[:, 1],
        "Euler3": e[:, 2],
        "MAD":    0.0,
        "BC":     200,
        "BS":     255,
    })
    return df


# ---------------------------------------------------------------------------
# Export helpers
# ---------------------------------------------------------------------------

def export_ctf(data: pd.DataFrame, filepath: str) -> None:
    """Write EBSD data back to a minimal CTF-format text file."""
    cols = ["Phase", "X", "Y", "Bands", "Error", "Euler1", "Euler2", "Euler3",
            "MAD", "BC", "BS"]
    for c in ["Bands", "Error"]:
        if c not in data.columns:
            data = data.copy()
            data[c] = 0
    out = data.reindex(columns=cols, fill_value=0)
    with open(filepath, "w") as fh:
        fh.write("Channel Text File\n")
        fh.write("Prj\tSAnTex_export\n")
        fh.write("Author\tSAnTex\n")
        fh.write("JobMode\tGrid\n")
        fh.write(f"XCells\t{int(data['X'].nunique())}\n")
        fh.write(f"YCells\t{int(data['Y'].nunique())}\n")
        xstep = float(np.diff(np.sort(data["X"].unique())).min()) if data["X"].nunique() > 1 else 1.0
        ystep = float(np.diff(np.sort(data["Y"].unique())).min()) if data["Y"].nunique() > 1 else 1.0
        fh.write(f"XStep\t{xstep:.4f}\n")
        fh.write(f"YStep\t{ystep:.4f}\n")
        fh.write("AcqE1\t0\nAcqE2\t0\nAcqE3\t0\n")
        fh.write("Euler angles refer to Sample Coordinate system (CS0)!\tMag\t0\t"
                 "Coverage\t0\tDevice\t0\tKV\t0\tTiltAngle\t0\tTiltAxis\t0\n")
        fh.write("Phases\t0\n")
        fh.write("\t".join(cols) + "\n")
        out.to_csv(fh, sep="\t", index=False, header=False, float_format="%.4f")


def export_ang(data: pd.DataFrame, filepath: str) -> None:
    """Write EBSD data in TSL .ang format (OIM-compatible)."""
    with open(filepath, "w") as fh:
        fh.write("# TEM_PIXperUM          1.000000\n")
        fh.write("# x-star                0.000000\n")
        fh.write("# y-star                0.000000\n")
        fh.write("# z-star                0.000000\n")
        fh.write("# WorkingDistance       0.000000\n")
        fh.write("#\n# OPERATOR:\tSAnTex\n#\n")
        fh.write("# PHASE 1\n#\n")
        fh.write("# GRID: SqrGrid\n")
        xstep = float(np.diff(np.sort(data["X"].unique())).min()) if data["X"].nunique() > 1 else 1.0
        ystep = float(np.diff(np.sort(data["Y"].unique())).min()) if data["Y"].nunique() > 1 else 1.0
        fh.write(f"# XSTEP: {xstep:.4f}\n")
        fh.write(f"# YSTEP: {ystep:.4f}\n")
        fh.write(f"# NCOLS_ODD: {data['X'].nunique()}\n")
        fh.write(f"# NCOLS_EVEN: {data['X'].nunique()}\n")
        fh.write(f"# NROWS: {data['Y'].nunique()}\n")
        fh.write("# HEADER: Done\n")
        # .ang columns: phi1 Phi phi2 x y IQ CI Phase SEM_signal Fit
        for _, row in data.iterrows():
            e1 = np.radians(row.get("Euler1", 0))
            e  = np.radians(row.get("Euler2", 0))
            e2 = np.radians(row.get("Euler3", 0))
            xi = row.get("X", 0)
            yi = row.get("Y", 0)
            bc = row.get("BC", 0)
            ci = row.get("MAD", 0)
            ph = int(row.get("Phase", 1))
            bs = row.get("BS", 0)
            fh.write(f"  {e1:.5f}  {e:.5f}  {e2:.5f}  {xi:.5f}  {yi:.5f}"
                     f"   {bc:.1f}   {ci:.4f}   {ph}   {bs:.1f}   0.0000\n")


def export_csv(data: pd.DataFrame, filepath: str) -> None:
    """Write EBSD data to a plain CSV."""
    data.to_csv(filepath, index=False, float_format="%.4f")


def export_hdf5(data: pd.DataFrame, filepath: str) -> None:
    """Write EBSD data to HDF5 (requires h5py)."""
    try:
        import h5py
    except ImportError:
        raise ImportError("h5py is required for HDF5 export: pip install h5py")
    with h5py.File(filepath, "w") as f:
        grp = f.create_group("EBSD/Data")
        for col in data.columns:
            try:
                grp.create_dataset(col, data=data[col].to_numpy())
            except Exception:
                pass


# ---------------------------------------------------------------------------
# Private rotation utilities (avoid importing santex.tensor inside tight loops)
# ---------------------------------------------------------------------------

def _euler_to_rotation_batch(phi1: np.ndarray, PHI: np.ndarray,
                               phi2: np.ndarray) -> np.ndarray:
    """Vectorised Bunge ZXZ → rotation matrix. Returns (N,3,3)."""
    c1, s1 = np.cos(phi1), np.sin(phi1)
    c,  s  = np.cos(PHI),  np.sin(PHI)
    c2, s2 = np.cos(phi2), np.sin(phi2)
    N = len(phi1)
    R = np.zeros((N, 3, 3))
    R[:, 0, 0] =  c1*c2 - s1*s2*c
    R[:, 0, 1] = -c1*s2 - s1*c2*c
    R[:, 0, 2] =  s1*s
    R[:, 1, 0] =  s1*c2 + c1*s2*c
    R[:, 1, 1] = -s1*s2 + c1*c2*c
    R[:, 1, 2] = -c1*s
    R[:, 2, 0] =  s2*s
    R[:, 2, 1] =  c2*s
    R[:, 2, 2] =  c
    return R


def _rotation_to_quaternion_batch(R: np.ndarray) -> np.ndarray:
    """Convert (N,3,3) rotation matrices to (N,4) unit quaternions [w,x,y,z]."""
    N = len(R)
    Q = np.zeros((N, 4))
    tr = R[:, 0, 0] + R[:, 1, 1] + R[:, 2, 2]
    Q[:, 0] = np.sqrt(np.maximum(0, (1 + tr) / 4))
    Q[:, 1] = np.sqrt(np.maximum(0, (1 + R[:, 0, 0] - R[:, 1, 1] - R[:, 2, 2]) / 4))
    Q[:, 2] = np.sqrt(np.maximum(0, (1 - R[:, 0, 0] + R[:, 1, 1] - R[:, 2, 2]) / 4))
    Q[:, 3] = np.sqrt(np.maximum(0, (1 - R[:, 0, 0] - R[:, 1, 1] + R[:, 2, 2]) / 4))
    Q[:, 1] = np.copysign(Q[:, 1], R[:, 2, 1] - R[:, 1, 2])
    Q[:, 2] = np.copysign(Q[:, 2], R[:, 0, 2] - R[:, 2, 0])
    Q[:, 3] = np.copysign(Q[:, 3], R[:, 1, 0] - R[:, 0, 1])
    norms = np.linalg.norm(Q, axis=1, keepdims=True)
    Q /= np.maximum(norms, 1e-12)
    return Q


def _quaternion_to_rotation(q: np.ndarray) -> np.ndarray:
    """Convert a single unit quaternion [w,x,y,z] to a (3,3) rotation matrix."""
    w, x, y, z = q
    return np.array([
        [1 - 2*(y*y + z*z),   2*(x*y - z*w),     2*(x*z + y*w)],
        [2*(x*y + z*w),       1 - 2*(x*x + z*z),  2*(y*z - x*w)],
        [2*(x*z - y*w),       2*(y*z + x*w),      1 - 2*(x*x + y*y)],
    ])


def _rotation_to_euler(R: np.ndarray) -> tuple[float, float, float]:
    """3×3 rotation matrix → Bunge ZXZ (φ₁, Φ, φ₂) in radians."""
    Phi = np.arccos(np.clip(R[2, 2], -1.0, 1.0))
    if np.abs(np.sin(Phi)) < 1e-10:
        phi1 = np.arctan2(R[0, 1], R[0, 0]) if Phi < 1.0 else np.arctan2(R[0, 1], -R[0, 0])
        phi2 = 0.0
    else:
        phi1 = np.arctan2(R[2, 0], -R[2, 1])
        phi2 = np.arctan2(R[0, 2],  R[1, 2])
    return phi1 % (2*np.pi), Phi, phi2 % (2*np.pi)


def _estimate_step(x: np.ndarray, y: np.ndarray) -> float:
    """Estimate the scan step size from the median NN distance of 200 random points."""
    from scipy.spatial import cKDTree
    n = min(len(x), 200)
    rng = np.random.default_rng(0)
    idx = rng.choice(len(x), n, replace=False)
    pts = np.column_stack([x[idx], y[idx]])
    tree = cKDTree(pts)
    dists, _ = tree.query(pts, k=2)          # k=2: self + nearest
    return float(np.median(dists[:, 1]))
