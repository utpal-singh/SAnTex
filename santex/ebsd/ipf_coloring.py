"""
Internal IPF colouring — pure NumPy, no orix dependency.

Architecture
------------
1. Comprehensive alias map: every HM / Schoenflies / short name → one of 11
   canonical Laue-group keys  (C1, Ci, C2h, D2h, C4h, D4h, S6, D3d,
   C6h, D6h, Th, Oh).
2. Elementary rotation-matrix generators (C2, C3, C4, C6, C3_111, …).
3. Group-closure algorithm — generates the full proper-rotation subgroup
   from a minimal set of generators.  Results are cached in a module-level
   dict so they are only ever built once per session.
4. Per-Laue-group fundamental-sector definition:
     • corner poles (unit vectors on the sphere)
     • corner RGB colours  (TSL / MTEX convention)
     • internal reference direction used for folding
5. Batch direction-folding: apply all symmetry ops + inversions, pick
   the equivalent direction closest to the reference.
6. Softmax IPF colouring: barycentric-like weighting by arccos distance
   from each corner pole  (T = 50, matching MTEX default sharpness).
7. Public :func:`ipf_map_colors` — same API as the former orix-based
   version in ebsd_operations.py.

References
----------
* MTEX source: ``mtex/geometry/symmetry/``
  https://github.com/mtex-toolbox/mtex
* TSL color-key convention: TSL OIM Analysis User Manual
* Softmax IPF coloring: Nolze & Hielscher (2016) J. Appl. Cryst.
"""

from __future__ import annotations
import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Canonical Laue-group keys and alias map
# ─────────────────────────────────────────────────────────────────────────────

#: The 11 canonical Laue-group keys used throughout this module.
LAUE_GROUPS: tuple[str, ...] = (
    "C1", "Ci",         # triclinic   (1,  -1)
    "C2h",              # monoclinic  (2/m)
    "D2h",              # orthorhombic (mmm)
    "C4h", "D4h",       # tetragonal  (4/m, 4/mmm)
    "S6",  "D3d",       # trigonal    (-3,  -3m)
    "C6h", "D6h",       # hexagonal   (6/m, 6/mmm)
    "Th",  "Oh",        # cubic       (m-3, m-3m)
)

_ALIAS: dict[str, str] = {
    # ── triclinic ──────────────────────────────────────────────────────────
    "1":    "C1",  "C1":  "C1",
    "-1":   "Ci",  "Ci":  "Ci",
    # ── monoclinic ─────────────────────────────────────────────────────────
    "2":    "C2h", "m":   "C2h", "2/m":  "C2h",
    "C2":   "C2h", "Cs":  "C2h", "C2h":  "C2h",
    # ── orthorhombic ───────────────────────────────────────────────────────
    "222":  "D2h", "mm2": "D2h", "mmm":  "D2h",
    "D2":   "D2h", "C2v": "D2h", "D2h":  "D2h",
    # ── tetragonal ─────────────────────────────────────────────────────────
    "4":    "C4h", "-4":  "C4h", "4/m":  "C4h",
    "C4":   "C4h", "S4":  "C4h", "C4h":  "C4h",
    "422":  "D4h", "4mm": "D4h", "-42m": "D4h", "-4m2": "D4h",
    "4/mmm":"D4h",
    "D4":   "D4h", "C4v": "D4h", "D2d":  "D4h", "D4h":  "D4h",
    # ── trigonal ───────────────────────────────────────────────────────────
    "3":    "S6",  "-3":  "S6",
    "C3":   "S6",  "C3i": "S6",  "S6":   "S6",
    "32":   "D3d", "3m":  "D3d", "-3m":  "D3d",
    "321":  "D3d", "3m1": "D3d", "-3m1": "D3d",
    "312":  "D3d", "31m": "D3d", "-31m": "D3d",
    "D3":   "D3d", "C3v": "D3d", "D3d":  "D3d",
    # ── hexagonal ──────────────────────────────────────────────────────────
    "6":    "C6h", "-6":  "C6h", "6/m":  "C6h",
    "C6":   "C6h", "C3h": "C6h", "C6h":  "C6h",
    "622":  "D6h", "6mm": "D6h", "-62m": "D6h", "-6m2": "D6h",
    "6/mmm":"D6h",
    "D6":   "D6h", "C6v": "D6h", "D3h":  "D6h", "D6h":  "D6h",
    # ── cubic ──────────────────────────────────────────────────────────────
    "23":   "Th",  "m-3": "Th",
    "T":    "Th",  "Th":  "Th",
    "432":  "Oh",  "-43m":"Oh",  "m-3m": "Oh",
    "O":    "Oh",  "Td":  "Oh",  "Oh":   "Oh",
}


def laue_group(name: str) -> str:
    """Return the canonical Laue-group key for *name*.

    Accepts Schoenflies (``'Oh'``, ``'D2h'``, …), Hermann–Mauguin short
    (``'m-3m'``, ``'mmm'``, …), and common variants (``'C3i'``, ``'-3m1'``,
    …).  Falls back to ``'D2h'`` (orthorhombic) for unknown input.
    """
    return _ALIAS.get(name, _ALIAS.get(name.strip(), "D2h"))


# ─────────────────────────────────────────────────────────────────────────────
# 2.  Elementary rotation-matrix generators
# ─────────────────────────────────────────────────────────────────────────────

_E = np.eye(3, dtype=np.float64)


def _Rz(deg: float) -> np.ndarray:
    a = np.deg2rad(deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]], dtype=np.float64)


def _C2n(n: np.ndarray) -> np.ndarray:
    """180° rotation about unit vector n:  R = 2 n⊗n − I."""
    nv = np.asarray(n, dtype=np.float64)
    nv = nv / np.linalg.norm(nv)
    return 2.0 * np.outer(nv, nv) - _E


# Named generators (proper rotations only)
_C2z = np.diag([-1., -1.,  1.]).astype(np.float64)
_C2y = np.diag([-1.,  1., -1.]).astype(np.float64)
_C2x = np.diag([ 1., -1., -1.]).astype(np.float64)

_C4z = np.array([[0., -1., 0.],
                  [1.,  0., 0.],
                  [0.,  0., 1.]], dtype=np.float64)

_C3z  = _Rz(120.)
_C3zm = _Rz(240.)
_C6z  = _Rz( 60.)

# 3-fold rotation about [111]: cycles x→y→z→x
#   R111 @ [1,0,0] = [0,1,0],  R111 @ [0,1,0] = [0,0,1],  R111 @ [0,0,1] = [1,0,0]
_C3_111 = np.array([[0., 0., 1.],
                     [1., 0., 0.],
                     [0., 1., 0.]], dtype=np.float64)


# ─────────────────────────────────────────────────────────────────────────────
# 3.  Group-closure algorithm
# ─────────────────────────────────────────────────────────────────────────────

def _group_closure(
    generators: list[np.ndarray],
    tol: float = 1e-9,
    max_ops: int = 96,
) -> np.ndarray:
    """Generate the complete group from *generators* by successive multiplication.

    Returns an (N, 3, 3) array of all unique proper-rotation matrices.
    Guaranteed to terminate: crystal point-group orders are ≤ 48 (Oh).
    """
    ops: list[np.ndarray] = [_E.copy()]
    for g in generators:
        g_arr = np.asarray(g, dtype=np.float64)
        if not any(np.allclose(g_arr, o, atol=tol) for o in ops):
            ops.append(g_arr)

    prev_len = 0
    while len(ops) > prev_len:
        prev_len = len(ops)
        snapshot = list(ops)          # snapshot so loop length is stable
        for a in snapshot:
            for b in snapshot:
                p = a @ b
                if not any(np.allclose(p, o, atol=tol) for o in ops):
                    ops.append(p)
                    if len(ops) >= max_ops:
                        return np.stack(ops, axis=0)

    return np.stack(ops, axis=0)      # (N, 3, 3)


# Module-level cache: str key → (N, 3, 3) ndarray
_SYM_CACHE: dict[str, np.ndarray] = {}


def get_sym_ops(name: str) -> np.ndarray:
    """Return (N, 3, 3) proper-rotation operators for the Laue group *name*.

    The operators are the *proper* rotations of the point group whose Laue
    class is *name*.  For direction-folding we also apply their negatives
    (inversions), which is handled in :func:`fold_directions`.

    Results are cached after the first call.
    """
    key = laue_group(name)
    if key in _SYM_CACHE:
        return _SYM_CACHE[key]

    if key == "C1":
        ops = _E[np.newaxis].copy()
    elif key == "Ci":
        ops = _E[np.newaxis].copy()
    elif key == "C2h":
        # monoclinic: {E, C2z}
        ops = _group_closure([_C2z])
    elif key == "D2h":
        # orthorhombic: {E, C2z, C2x, C2y}
        ops = _group_closure([_C2z, _C2x])
    elif key == "C4h":
        # tetragonal 4/m: {E, C4z, C2z, C4zm}
        ops = _group_closure([_C4z])
    elif key == "D4h":
        # tetragonal 4/mmm: 8-element dihedral group D4
        ops = _group_closure([_C4z, _C2x])
    elif key == "S6":
        # trigonal -3 (C3i): {E, C3z, C3zm}
        ops = _group_closure([_C3z])
    elif key == "D3d":
        # trigonal -3m: 6-element D3 — {E, C3z, C3zm, C2x, C2@120, C2@240}
        ops = _group_closure([_C3z, _C2x])
    elif key == "C6h":
        # hexagonal 6/m: {E, C6z, C3z, C2z, C3zm, C6zm}
        ops = _group_closure([_C6z])
    elif key == "D6h":
        # hexagonal 6/mmm: 12-element D6
        ops = _group_closure([_C6z, _C2x])
    elif key == "Th":
        # cubic m-3: 12-element chiral-tetrahedral T
        ops = _group_closure([_C2z, _C2x, _C3_111])
    elif key == "Oh":
        # cubic m-3m: 24-element chiral-octahedral O
        ops = _group_closure([_C4z, _C3_111])
    else:
        ops = _E[np.newaxis].copy()

    _SYM_CACHE[key] = ops
    return ops


# ─────────────────────────────────────────────────────────────────────────────
# 4.  Fundamental-sector data (corner poles + colours + reference direction)
# ─────────────────────────────────────────────────────────────────────────────

_s2 = np.sqrt(2.)
_s3 = np.sqrt(3.)

# For each Laue group:
#   "corners"  – (K, 3) normalised corner unit vectors on the sphere
#   "colors"   – (K, 3) RGB corner colours  (TSL / MTEX convention)
#   "ref"      – (3,)   reference direction *inside* the sector; used by
#                fold_directions() to pick the canonical equivalent direction.
#
# TSL / MTEX colour convention:
#   [001] / c-axis  → blue  (0, 0, 1)
#   [100] / a1-axis → red   (1, 0, 0)  [or green for some groups]
#   [110] / edge    → green (0, 1, 0)
#   [111] / body    → red   (1, 0, 0)  [cubic only]
#
_SECTOR: dict[str, dict] = {

    # ── triclinic ───────────────────────────────────────────────────────────
    "C1": {
        "corners": np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]]),
        "colors":  np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]]),
        "ref":     np.array([1., 1., 1.]) / _s3,
    },
    "Ci": {
        "corners": np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]]),
        "colors":  np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]]),
        "ref":     np.array([1., 1., 1.]) / _s3,
    },

    # ── monoclinic ──────────────────────────────────────────────────────────
    # Unique axis is z; fundamental sector is upper half-hemisphere (180°).
    # Two corners suffice for a reasonable colour gradient.
    "C2h": {
        "corners": np.array([[0., 0., 1.], [1., 0., 0.]]),
        "colors":  np.array([[0., 0., 1.], [1., 0., 0.]]),
        "ref":     np.array([0., 0., 1.]),
    },

    # ── orthorhombic ────────────────────────────────────────────────────────
    # IPF triangle: [001] – [010] – [100]  (upper octant, z≥y≥0, x≥0)
    # TSL: [001]=blue, [010]=green, [100]=red
    "D2h": {
        "corners": np.array([[0., 0., 1.], [0., 1., 0.], [1., 0., 0.]]),
        "colors":  np.array([[0., 0., 1.], [0., 1., 0.], [1., 0., 0.]]),
        "ref":     np.array([1., 1., 1.]) / _s3,
    },

    # ── tetragonal ──────────────────────────────────────────────────────────
    "C4h": {
        "corners": np.array([[0., 0., 1.], [1., 0., 0.]]),
        "colors":  np.array([[0., 0., 1.], [1., 0., 0.]]),
        "ref":     np.array([0., 0., 1.]),
    },
    # IPF triangle: [001] – [100] – [110]/√2
    # TSL: [001]=blue, [100]=red, [110]=green
    "D4h": {
        "corners": np.array([
            [0., 0., 1.],
            [1., 0., 0.],
            [1./_s2, 1./_s2, 0.],
        ]),
        "colors": np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]]),
        "ref":    np.array([1., 0.5, 2.]) / np.linalg.norm([1., 0.5, 2.]),
    },

    # ── trigonal ────────────────────────────────────────────────────────────
    "S6": {
        "corners": np.array([[0., 0., 1.], [1., 0., 0.]]),
        "colors":  np.array([[0., 0., 1.], [1., 0., 0.]]),
        "ref":     np.array([0., 0., 1.]),
    },
    # IPF triangle: [001] – [100] – [0.5, √3/2, 0]  (azimuth 0° → 60°)
    # TSL: [001]=blue, [100]=red, [2̄110]=green
    "D3d": {
        "corners": np.array([
            [0.,   0.,        1.],
            [1.,   0.,        0.],
            [0.5,  _s3 / 2.,  0.],      # 60° from x in basal plane
        ]),
        "colors": np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]]),
        "ref":    np.array([np.cos(np.deg2rad(30.)),
                             np.sin(np.deg2rad(30.)), 1.]),
    },

    # ── hexagonal ───────────────────────────────────────────────────────────
    "C6h": {
        "corners": np.array([[0., 0., 1.], [1., 0., 0.]]),
        "colors":  np.array([[0., 0., 1.], [1., 0., 0.]]),
        "ref":     np.array([0., 0., 1.]),
    },
    # IPF triangle: [001] – [100] – [√3/2, 1/2, 0]  (azimuth 0° → 30°)
    # TSL: [001]=blue, [100]=red, [11̄20]=green
    "D6h": {
        "corners": np.array([
            [0.,       0.,   1.],
            [1.,       0.,   0.],
            [_s3/2.,   0.5,  0.],       # 30° from x in basal plane
        ]),
        "colors": np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]]),
        "ref":    np.array([np.cos(np.deg2rad(15.)),
                             np.sin(np.deg2rad(15.)), 1.]),
    },

    # ── cubic ───────────────────────────────────────────────────────────────
    # Th (m-3): IPF triangle similar to Oh; 12 proper rotations (T group)
    "Th": {
        "corners": np.array([
            [0.,      0.,      1.],
            [1./_s2,  1./_s2,  0.],     # [110]
            [1./_s3,  1./_s3,  1./_s3], # [111]
        ]),
        "colors": np.array([[0., 0., 1.], [0., 1., 0.], [1., 0., 0.]]),
        "ref":    np.array([1., 2., 3.]) / np.linalg.norm([1., 2., 3.]),
    },
    # Oh (m-3m): standard cubic IPF triangle
    # [001]=blue, [011]/√2=green, [111]/√3=red
    "Oh": {
        "corners": np.array([
            [0.,      0.,      1.],       # [001] → blue
            [0.,      1./_s2,  1./_s2],   # [011] → green
            [1./_s3,  1./_s3,  1./_s3],  # [111] → red
        ]),
        "colors": np.array([[0., 0., 1.], [0., 1., 0.], [1., 0., 0.]]),
        "ref":    np.array([1., 2., 3.]) / np.linalg.norm([1., 2., 3.]),
    },
}

# Normalise corner vectors and reference directions in-place
for _key, _sec in _SECTOR.items():
    _c = _sec["corners"].astype(np.float64)
    _norms = np.linalg.norm(_c, axis=1, keepdims=True)
    _sec["corners"] = np.where(_norms > 0, _c / _norms, _c)
    _r = np.asarray(_sec["ref"], dtype=np.float64)
    _sec["ref"] = _r / np.linalg.norm(_r)
del _key, _sec, _c, _norms, _r


# ─────────────────────────────────────────────────────────────────────────────
# 5.  Batch direction folding
# ─────────────────────────────────────────────────────────────────────────────

def fold_directions(V: np.ndarray, sym_key: str) -> np.ndarray:
    """Map each direction in *V* to the canonical fundamental-sector representative.

    Parameters
    ----------
    V       : (N, 3) unit vectors (crystal frame)
    sym_key : canonical Laue-group key from :func:`laue_group`

    Returns
    -------
    (N, 3) unit vectors whose azimuth / polar angles lie in the fundamental
    sector of *sym_key*.

    Algorithm
    ---------
    For each pixel:
    1. Apply all M proper-rotation operators → M equivalent directions.
    2. Append their negatives (because crystal directions are axial, v ≡ −v).
    3. Pick the equivalent direction whose dot-product with the internal
       reference direction is maximal.  The reference lives well inside the
       fundamental sector, so the dot-product maximum uniquely selects the
       canonical representative.
    """
    ops = get_sym_ops(sym_key)                                # (M, 3, 3)
    ref = _SECTOR.get(sym_key, _SECTOR["D2h"])["ref"]        # (3,)

    # equiv[i, j, :] = ops[j] @ V[i]   shape: (N, M, 3)
    equiv = np.einsum("jkl,il->ijk", ops, V)

    # Axial directions: include negatives — (N, 2M, 3)
    equiv = np.concatenate([equiv, -equiv], axis=1)

    # Score: projection onto reference direction — (N, 2M)
    scores = np.einsum("ijk,k->ij", equiv, ref)

    best = np.argmax(scores, axis=1)                          # (N,)
    N = len(V)
    canonical = equiv[np.arange(N), best]                    # (N, 3)

    # Ensure unit length (numerical safety)
    norms = np.linalg.norm(canonical, axis=1, keepdims=True)
    canonical = np.where(norms > 1e-15, canonical / norms, canonical)
    return canonical


# ─────────────────────────────────────────────────────────────────────────────
# 6.  Softmax IPF colour
# ─────────────────────────────────────────────────────────────────────────────

def _ipf_color_from_sector(
    v_canonical: np.ndarray,
    sector_key: str,
    T: float = 50.0,
) -> np.ndarray:
    """Assign per-pixel RGB via softmax distance weighting from corner poles.

    Parameters
    ----------
    v_canonical : (N, 3) unit vectors in the fundamental sector
    sector_key  : canonical Laue-group key
    T           : softmax temperature — larger → sharper colour transitions
                  (50 matches MTEX default)

    Returns
    -------
    (N, 3) RGB in [0, 1]
    """
    sec = _SECTOR.get(sector_key, _SECTOR["D2h"])
    corners = sec["corners"]    # (K, 3)
    colors  = sec["colors"]     # (K, 3)

    # Angle between each direction and each corner: (N, K)
    cos_ang = np.clip(v_canonical @ corners.T, -1.0, 1.0)
    angles  = np.arccos(cos_ang)                             # radians

    # Softmax weights:  w_k = exp(−T · θ_k²)
    weights = np.exp(-T * angles ** 2)                       # (N, K)
    weights /= weights.sum(axis=1, keepdims=True) + 1e-30    # normalise rows

    # Weighted sum of corner colours
    rgb = weights @ colors                                   # (N, 3)
    return np.clip(rgb, 0.0, 1.0)


# ─────────────────────────────────────────────────────────────────────────────
# 7.  Private rotation utilities (self-contained — no cross-module import)
# ─────────────────────────────────────────────────────────────────────────────

def _euler_to_rotation_batch(
    phi1: np.ndarray, PHI: np.ndarray, phi2: np.ndarray
) -> np.ndarray:
    """Vectorised Bunge ZXZ → rotation matrix.

    Parameters: three (N,) arrays of angles in **radians**.
    Returns: (N, 3, 3) rotation matrices.
    """
    c1, s1 = np.cos(phi1), np.sin(phi1)
    c,  s  = np.cos(PHI),  np.sin(PHI)
    c2, s2 = np.cos(phi2), np.sin(phi2)
    N = len(phi1)
    R = np.zeros((N, 3, 3), dtype=np.float64)
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


# ─────────────────────────────────────────────────────────────────────────────
# 8.  Public API — ipf_map_colors
# ─────────────────────────────────────────────────────────────────────────────

def ipf_map_colors(
    data,                               # pd.DataFrame
    phase_index: int = 1,
    direction: str = "ND",
    crystal_symmetry: str = "D2h",
) -> dict | None:
    """Compute per-pixel IPF colour for *phase_index* — no orix dependency.

    Parameters
    ----------
    data            : EBSD DataFrame with columns Phase, X, Y,
                      Euler1, Euler2, Euler3 (angles in degrees).
    phase_index     : phase number to colour (pixels with Phase == 0 are
                      excluded).
    direction       : sample reference direction —
                      ``'ND'`` (normal, Z), ``'RD'`` (rolling, X),
                      ``'TD'`` (transverse, Y).
    crystal_symmetry: any accepted symmetry name — Schoenflies (``'Oh'``,
                      ``'D2h'``, …), HM (``'m-3m'``, ``'mmm'``, …), or
                      common abbreviations (``'C3i'``, ``'-3m1'``, …).

    Returns
    -------
    dict with keys:

    ============  ========================================
    ``x``, ``y``  pixel coordinates
    ``r``, ``g``, ``b``  float RGB in [0, 1]
    ``rgb_hex``   list of ``'#rrggbb'`` strings
    ============  ========================================

    Returns ``None`` if *phase_index* has no indexed pixels.
    """
    _DIR = {
        "ND": np.array([0., 0., 1.]),
        "RD": np.array([1., 0., 0.]),
        "TD": np.array([0., 1., 0.]),
    }
    sample_dir = _DIR.get(direction.upper(), _DIR["ND"])

    phase_data = data[data["Phase"] == phase_index]
    if len(phase_data) == 0:
        return None

    # Euler angles → rotation matrices (N, 3, 3)
    phi1 = np.deg2rad(phase_data["Euler1"].to_numpy(dtype=float))
    PHI  = np.deg2rad(phase_data["Euler2"].to_numpy(dtype=float))
    phi2 = np.deg2rad(phase_data["Euler3"].to_numpy(dtype=float))
    R = _euler_to_rotation_batch(phi1, PHI, phi2)    # (N, 3, 3)

    # Crystal direction corresponding to sample_dir:
    #   d_crystal = R^T @ d_sample  (passive rotation convention)
    crystal_dirs = np.einsum("nij,j->ni", R, sample_dir)    # (N, 3)

    # Normalise (handles rounding errors at poles)
    norms = np.linalg.norm(crystal_dirs, axis=1, keepdims=True)
    crystal_dirs = np.where(norms > 1e-15, crystal_dirs / norms, crystal_dirs)

    # Map to fundamental sector
    sym_key = laue_group(crystal_symmetry)
    canonical = fold_directions(crystal_dirs, sym_key)       # (N, 3)

    # Assign RGB
    rgb = _ipf_color_from_sector(canonical, sym_key)         # (N, 3)

    hex_colors = [
        "#{:02x}{:02x}{:02x}".format(
            int(r * 255), int(g * 255), int(b * 255))
        for r, g, b in rgb
    ]

    return {
        "x":       phase_data["X"].to_numpy(),
        "y":       phase_data["Y"].to_numpy(),
        "r":       rgb[:, 0],
        "g":       rgb[:, 1],
        "b":       rgb[:, 2],
        "rgb_hex": hex_colors,
    }
