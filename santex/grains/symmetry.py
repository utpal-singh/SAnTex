"""
Crystal symmetry operations for misorientation calculations.

Covers all 32 crystallographic point groups organised by Laue group.
Each Laue group is represented by the proper rotation matrices (det = +1)
of its maximal proper rotation subgroup.

MTEX-compatible aliases are supported (case-insensitive).  The 11 Laue
groups and their proper-rotation orders are:

  Laue   Schoenflies  System            n_ops
  -----  -----------  ----------------  -----
  -1     Ci           triclinic             1
  2/m    C2h          monoclinic            2
  mmm    D2h          orthorhombic          4
  4/m    C4h          tetragonal (low)      4
  4/mmm  D4h          tetragonal (high)     8
  -3     S6           trigonal (low)        3
  -3m    D3d          trigonal (high)       6
  6/m    C6h          hexagonal (low)       6
  6/mmm  D6h          hexagonal (high)     12
  m-3    Th           cubic (low)          12
  m-3m   Oh           cubic (high)         24

References
----------
Hielscher & Schaeben (2008), J. Appl. Cryst. 41, 1024-1037
MTEX toolbox: https://mtex-toolbox.github.io
International Tables for Crystallography Vol. A
"""

from __future__ import annotations
import itertools
import numpy as np

# ---------------------------------------------------------------------------
# Low-level rotation matrix helpers
# ---------------------------------------------------------------------------

def _Rz(a: float) -> np.ndarray:
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]], dtype=np.float64)


def _Ry(a: float) -> np.ndarray:
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, 0., s], [0., 1., 0.], [-s, 0., c]], dtype=np.float64)


def _Rx(a: float) -> np.ndarray:
    c, s = np.cos(a), np.sin(a)
    return np.array([[1., 0., 0.], [0., c, -s], [0., s, c]], dtype=np.float64)


# ---------------------------------------------------------------------------
# Proper rotation group generators  (one function per Laue proper subgroup)
# ---------------------------------------------------------------------------

def _C1_sym() -> np.ndarray:
    """Triclinic C1 — 1 proper rotation (identity)."""
    return np.eye(3, dtype=np.float64)[np.newaxis]


def _C2_sym() -> np.ndarray:
    """Monoclinic C2 — 2 proper rotations: E, C2y (b-unique)."""
    return np.array([
        np.eye(3),
        np.diag([-1., 1., -1.]),   # C2 about y
    ], dtype=np.float64)


def _C3_sym() -> np.ndarray:
    """Trigonal C3 — 3 proper rotations: E, C3z, C3z²."""
    return np.array([_Rz(k * 2.0 * np.pi / 3.0) for k in range(3)], dtype=np.float64)


def _C4_sym() -> np.ndarray:
    """Tetragonal C4 — 4 proper rotations: E, C4z, C2z, C4z³."""
    return np.array([_Rz(k * np.pi / 2.0) for k in range(4)], dtype=np.float64)


def _C6_sym() -> np.ndarray:
    """Hexagonal C6 — 6 proper rotations: E, C6z, C3z, C2z, C3z², C6z⁵."""
    return np.array([_Rz(k * np.pi / 3.0) for k in range(6)], dtype=np.float64)


def _D2_sym() -> np.ndarray:
    """Orthorhombic D2 — 4 proper rotations: E, C2z, C2y, C2x."""
    return np.array([
        np.eye(3),
        np.diag([-1., -1.,  1.]),   # C2z
        np.diag([-1.,  1., -1.]),   # C2y
        np.diag([ 1., -1., -1.]),   # C2x
    ], dtype=np.float64)


def _D3_sym() -> np.ndarray:
    """Trigonal D3 — 6 proper rotations: 3 C3 + 3 C2 in-plane."""
    ops = []
    for k in range(3):
        a = k * 2.0 * np.pi / 3.0
        ops.append(_Rz(a))
        c2, s2 = np.cos(2.0 * a), np.sin(2.0 * a)
        ops.append(np.array([[c2, s2, 0.], [s2, -c2, 0.], [0., 0., -1.]], dtype=np.float64))
    return np.array(ops, dtype=np.float64)


def _D4_sym() -> np.ndarray:
    """Tetragonal D4 — 8 proper rotations: 4 C4 + 4 C2."""
    ops = [_Rz(k * np.pi / 2.0) for k in range(4)]
    for k in range(4):
        a = k * np.pi / 4.0
        c2, s2 = np.cos(2.0 * a), np.sin(2.0 * a)
        ops.append(np.array([[c2, s2, 0.], [s2, -c2, 0.], [0., 0., -1.]], dtype=np.float64))
    return np.array(ops, dtype=np.float64)


def _D6_sym() -> np.ndarray:
    """Hexagonal D6 — 12 proper rotations: 6 C6 + 6 C2."""
    ops = [_Rz(k * np.pi / 3.0) for k in range(6)]
    for k in range(6):
        a = k * np.pi / 3.0
        c2, s2 = np.cos(2.0 * a), np.sin(2.0 * a)
        ops.append(np.array([[c2, s2, 0.], [s2, -c2, 0.], [0., 0., -1.]], dtype=np.float64))
    return np.array(ops, dtype=np.float64)


def _proper_rots_from_perms_signs() -> np.ndarray:
    """All 48 permutation+sign matrices; return 24 with det = +1 (cubic O)."""
    ops = []
    for perm in itertools.permutations([0, 1, 2]):
        for signs in itertools.product([-1, 1], repeat=3):
            M = np.zeros((3, 3))
            for row, (col, s) in enumerate(zip(perm, signs)):
                M[row, col] = s
            if round(np.linalg.det(M)) == 1:
                ops.append(M)
    return np.array(ops, dtype=np.float64)


def _T_sym() -> np.ndarray:
    """
    Chiral tetrahedral T (23 / m-3) — 12 proper rotations.

    Elements: E, 3 C2 about face-normals [100]/[010]/[001],
              8 C3 about body diagonals [111]/[1-1-1]/[-11-1]/[-1-11].

    Selected from O (24 ops) by excluding the 6 C4 and 6 C2-edge operations.
    """
    O = _proper_rots_from_perms_signs()
    T_ops = []
    for M in O:
        tr = np.trace(M)
        if abs(tr - 3.0) < 0.1:     # identity   (trace = 3)
            T_ops.append(M)
        elif abs(tr) < 0.1:          # C3 / C3²   (trace = 0)
            T_ops.append(M)
        elif abs(tr + 1.0) < 0.1:    # C2         (trace = -1)
            # Distinguish face-axis C2 (in T) from edge-axis C2 (not in T)
            vals, vecs = np.linalg.eig(M)
            idx = int(np.argmin(np.abs(np.real(vals) - 1.0)))
            axis = np.abs(np.real(vecs[:, idx]))
            # Face axis [1,0,0]-type: second-largest component ≈ 0
            axis_sorted = np.sort(axis)[::-1]
            if axis_sorted[1] < 0.15:
                T_ops.append(M)
    return np.array(T_ops, dtype=np.float64)


def _O_sym() -> np.ndarray:
    """Cubic O (432 / m-3m) — 24 proper rotations."""
    return _proper_rots_from_perms_signs()


# ---------------------------------------------------------------------------
# Backward-compatible aliases (match old function names)
# ---------------------------------------------------------------------------

_CUBIC_SYM      = _O_sym()
_hex_symmetry   = _D6_sym
_tetra_symmetry = _D4_sym
_ortho_symmetry = _D2_sym
_trig_symmetry  = _D3_sym
_mono_symmetry  = _C2_sym
_tri_symmetry   = _C1_sym


# ---------------------------------------------------------------------------
# The 11 Laue-group proper-rotation subgroups
# ---------------------------------------------------------------------------

_LAUE_PROPER_SYM: dict[str, np.ndarray] = {
    "-1":    _C1_sym(),    # triclinic
    "2/m":   _C2_sym(),    # monoclinic
    "mmm":   _D2_sym(),    # orthorhombic
    "4/m":   _C4_sym(),    # tetragonal (low)
    "4/mmm": _D4_sym(),    # tetragonal (high)
    "-3":    _C3_sym(),    # trigonal (low)
    "-3m":   _D3_sym(),    # trigonal (high) / rhombohedral
    "6/m":   _C6_sym(),    # hexagonal (low)
    "6/mmm": _D6_sym(),    # hexagonal (high)
    "m-3":   _T_sym(),     # cubic (low)
    "m-3m":  _O_sym(),     # cubic (high)
}


# ---------------------------------------------------------------------------
# Full crystallographic symmetry metadata
# ---------------------------------------------------------------------------

CRYSTAL_SYMMETRY_INFO: dict[str, dict] = {
    "-1": {
        "display":        "Triclinic (-1)",
        "crystal_system": "triclinic",
        "laue_group":     "-1",
        "laue_id":        1,
        "proper_order":   1,
        "schoenflies":    "Ci",
        "hm_symbol":      "-1",
        "point_groups":   ["1 (C1)", "-1 (Ci / S2)"],
        "aliases":        ["1", "-1", "c1", "ci", "s2",
                           "triclinic", "anorthic"],
        "minerals":       ["plagioclase", "microcline", "kyanite", "axinite", "wollastonite"],
    },
    "2/m": {
        "display":        "Monoclinic (2/m)",
        "crystal_system": "monoclinic",
        "laue_group":     "2/m",
        "laue_id":        2,
        "proper_order":   2,
        "schoenflies":    "C2h",
        "hm_symbol":      "2/m",
        "point_groups":   ["2 (C2)", "m (Cs / C1h)", "2/m (C2h)"],
        "aliases":        ["2", "m", "2/m", "c2", "cs", "c1h", "c2h",
                           "monoclinic", "2/m11", "12/m1", "112/m"],
        "minerals":       ["gypsum", "orthoclase", "hornblende", "augite",
                           "muscovite", "epidote", "diopside"],
    },
    "mmm": {
        "display":        "Orthorhombic (mmm)",
        "crystal_system": "orthorhombic",
        "laue_group":     "mmm",
        "laue_id":        3,
        "proper_order":   4,
        "schoenflies":    "D2h",
        "hm_symbol":      "mmm",
        "point_groups":   ["222 (D2 / V)", "mm2 (C2v)", "mmm (D2h / Vh)"],
        "aliases":        ["222", "mm2", "2mm", "m2m", "mmm",
                           "d2", "v", "c2v", "d2h", "vh",
                           "orthorhombic"],
        "minerals":       ["olivine", "forsterite", "enstatite", "fayalite",
                           "topaz", "aragonite", "sulfur", "barite"],
    },
    "4/m": {
        "display":        "Tetragonal low (4/m)",
        "crystal_system": "tetragonal",
        "laue_group":     "4/m",
        "laue_id":        4,
        "proper_order":   4,
        "schoenflies":    "C4h",
        "hm_symbol":      "4/m",
        "point_groups":   ["4 (C4)", "-4 (S4)", "4/m (C4h)"],
        "aliases":        ["4", "-4", "4/m", "c4", "s4", "c4h",
                           "tetragonal-low"],
        "minerals":       ["wulfenite", "scheelite", "cahnite"],
    },
    "4/mmm": {
        "display":        "Tetragonal (4/mmm)",
        "crystal_system": "tetragonal",
        "laue_group":     "4/mmm",
        "laue_id":        5,
        "proper_order":   8,
        "schoenflies":    "D4h",
        "hm_symbol":      "4/mmm",
        "point_groups":   ["422 (D4)", "4mm (C4v)", "-42m (D2d / Vd)", "4/mmm (D4h)"],
        "aliases":        ["422", "4mm", "-42m", "-4m2", "4/mmm",
                           "d4", "c4v", "d2d", "vd", "d4h",
                           "tetragonal", "42"],
        "minerals":       ["zircon", "rutile", "cassiterite", "vesuvianite",
                           "apophyllite", "indium", "tin"],
    },
    "-3": {
        "display":        "Trigonal low (-3)",
        "crystal_system": "trigonal",
        "laue_group":     "-3",
        "laue_id":        6,
        "proper_order":   3,
        "schoenflies":    "S6",
        "hm_symbol":      "-3",
        "point_groups":   ["3 (C3)", "-3 (S6 / C3i)"],
        "aliases":        ["3", "-3", "c3", "s6", "c3i",
                           "trigonal-low"],
        "minerals":       ["phenakite", "dioptase", "jarosite"],
    },
    "-3m": {
        "display":        "Trigonal / Rhombohedral (-3m)",
        "crystal_system": "trigonal",
        "laue_group":     "-3m",
        "laue_id":        7,
        "proper_order":   6,
        "schoenflies":    "D3d",
        "hm_symbol":      "-3m",
        "point_groups":   ["32 (D3)", "3m (C3v)", "-3m1 (D3d)", "-31m (D3d)"],
        "aliases":        ["32", "3m", "-3m", "-3m1", "-31m",
                           "d3", "c3v", "d3d",
                           "trigonal", "rhombohedral"],
        "minerals":       ["calcite", "dolomite", "quartz", "tourmaline",
                           "corundum", "hematite", "bismuth", "antimony"],
    },
    "6/m": {
        "display":        "Hexagonal low (6/m)",
        "crystal_system": "hexagonal",
        "laue_group":     "6/m",
        "laue_id":        8,
        "proper_order":   6,
        "schoenflies":    "C6h",
        "hm_symbol":      "6/m",
        "point_groups":   ["6 (C6)", "-6 (C3h)", "6/m (C6h)"],
        "aliases":        ["6", "-6", "6/m", "c6", "c3h", "c6h",
                           "hexagonal-low"],
        "minerals":       ["nepheline (6/m)", "apatite (6/m)", "covellite"],
    },
    "6/mmm": {
        "display":        "Hexagonal (6/mmm)",
        "crystal_system": "hexagonal",
        "laue_group":     "6/mmm",
        "laue_id":        9,
        "proper_order":   12,
        "schoenflies":    "D6h",
        "hm_symbol":      "6/mmm",
        "point_groups":   ["622 (D6)", "6mm (C6v)", "-6m2 (D3h)", "6/mmm (D6h)"],
        "aliases":        ["622", "6mm", "-6m2", "-62m", "6/mmm",
                           "d6", "c6v", "d3h", "d6h",
                           "hexagonal", "62"],
        "minerals":       ["beryl", "molybdenite", "graphite", "ice Ih",
                           "apatite (6/mmm)", "titanium", "magnesium",
                           "zinc", "zirconium", "cobalt"],
    },
    "m-3": {
        "display":        "Cubic low (m-3)",
        "crystal_system": "cubic",
        "laue_group":     "m-3",
        "laue_id":        10,
        "proper_order":   12,
        "schoenflies":    "Th",
        "hm_symbol":      "m-3",
        "point_groups":   ["23 (T)", "m-3 (Th)"],
        "aliases":        ["23", "m-3", "m3", "t", "th",
                           "cubic-low", "pyrite"],
        "minerals":       ["pyrite", "sperrylite", "ullmannite",
                           "cobaltite", "skutterudite"],
    },
    "m-3m": {
        "display":        "Cubic (m-3m)",
        "crystal_system": "cubic",
        "laue_group":     "m-3m",
        "laue_id":        11,
        "proper_order":   24,
        "schoenflies":    "Oh",
        "hm_symbol":      "m-3m",
        "point_groups":   ["432 (O)", "-43m (Td)", "m-3m (Oh)"],
        "aliases":        ["432", "-43m", "43m", "m-3m", "m3m",
                           "o", "td", "oh",
                           "cubic", "isometric", "43"],
        "minerals":       ["iron", "copper", "aluminum", "nickel",
                           "garnet", "spinel", "magnetite",
                           "diamond", "halite", "galena", "fluorite"],
    },
}


# ---------------------------------------------------------------------------
# Point-group → Laue-group lookup  (all 32 point groups + MTEX aliases)
# ---------------------------------------------------------------------------

_PG_TO_LAUE: dict[str, str] = {}

# Populate from CRYSTAL_SYMMETRY_INFO: every alias → its Laue group key
for _laue_key, _info in CRYSTAL_SYMMETRY_INFO.items():
    for _alias in _info["aliases"]:
        _PG_TO_LAUE[_alias.lower()] = _laue_key

# Additional MTEX findsymmetry.m aliases not already covered
_EXTRA_ALIASES: dict[str, str] = {
    # Triclinic
    "c1":          "-1",
    # Monoclinic settings
    "2/m11":       "2/m",
    "12/m1":       "2/m",
    "112/m":       "2/m",
    # Orthorhombic
    "vh":          "mmm",
    "v":           "mmm",
    # Orthorhombic mineral names frequently used in geoscience
    "olivine":     "mmm",
    "forsterite":  "mmm",
    "enstatite":   "mmm",
    "fayalite":    "mmm",
    # Tetragonal: 4/m vs 4/mmm
    "4/m":         "4/m",
    "4/mmm":       "4/mmm",
    # Trigonal
    "-3m1":        "-3m",
    "-31m":        "-3m",
    "3m1":         "-3m",
    "31m":         "-3m",
    # Cubic low
    "m-3":         "m-3",
    "m3":          "m-3",
    # Cubic high
    "m-3m":        "m-3m",
    "m3m":         "m-3m",
    # General crystal-system words (legacy)
    "unknown":     "-1",
    "anorthic":    "-1",
}
for _k, _v in _EXTRA_ALIASES.items():
    _PG_TO_LAUE.setdefault(_k.lower(), _v)


# ---------------------------------------------------------------------------
# Main public SYMMETRY dict
# (all lowercase keys → proper rotation matrix array)
# ---------------------------------------------------------------------------

SYMMETRY: dict[str, np.ndarray] = {
    key: _LAUE_PROPER_SYM[laue]
    for key, laue in _PG_TO_LAUE.items()
}


# ---------------------------------------------------------------------------
# Symmetry-pair cache
# ---------------------------------------------------------------------------

_SYM_PAIRS: dict[str, np.ndarray] = {}


def _get_sym_pairs(key: str) -> np.ndarray:
    key = key.lower().strip()
    if key not in _SYM_PAIRS:
        _SYM_PAIRS[key] = SYMMETRY.get(key, _LAUE_PROPER_SYM["-1"])
    return _SYM_PAIRS[key]


# ---------------------------------------------------------------------------
# Public API — symmetry lookup
# ---------------------------------------------------------------------------

def get_symmetry(name: str) -> np.ndarray:
    """
    Return the proper rotation matrices for the named crystal symmetry.

    Accepts any of the 32 point-group HM symbols, Schoenflies symbols,
    Laue-group symbols, crystal-system words, or MTEX aliases
    (all case-insensitive).

    Parameters
    ----------
    name : str
        E.g. ``"m-3m"``, ``"cubic"``, ``"Oh"``, ``"D6h"``, ``"6/mmm"``,
        ``"-43m"``, ``"orthorhombic"``, ``"trigonal"``, …

    Returns
    -------
    ops : (N, 3, 3) ndarray  — proper rotation matrices
    """
    key = name.lower().strip()
    laue = _PG_TO_LAUE.get(key, "-1")
    return _LAUE_PROPER_SYM[laue]


def get_laue_group(name: str) -> str:
    """
    Return the Laue-group HM symbol for the given symmetry name.

    Returns ``"-1"`` (triclinic) for unrecognised names.
    """
    return _PG_TO_LAUE.get(name.lower().strip(), "-1")


def get_info(name: str) -> dict:
    """
    Return the full metadata dict for the Laue group that contains *name*.

    Keys: ``display``, ``crystal_system``, ``laue_group``, ``laue_id``,
    ``proper_order``, ``schoenflies``, ``hm_symbol``, ``point_groups``,
    ``aliases``, ``minerals``.
    """
    laue = get_laue_group(name)
    return CRYSTAL_SYMMETRY_INFO[laue]


def list_symmetries(verbose: bool = False) -> list[str]:
    """
    Return sorted list of all recognised symmetry aliases.

    Parameters
    ----------
    verbose : bool
        If True, return display strings (one per Laue group) instead of aliases.
    """
    if verbose:
        return [v["display"] for v in CRYSTAL_SYMMETRY_INFO.values()]
    return sorted(_PG_TO_LAUE.keys())


def list_laue_groups() -> list[dict]:
    """
    Return a list of dicts (one per Laue group) sorted by ``laue_id``.

    Each dict has keys: ``laue_group``, ``display``, ``crystal_system``,
    ``laue_id``, ``proper_order``, ``schoenflies``.
    """
    return sorted(
        [
            {k: v[k] for k in ("laue_group", "display", "crystal_system",
                                "laue_id", "proper_order", "schoenflies")}
            for v in CRYSTAL_SYMMETRY_INFO.values()
        ],
        key=lambda d: d["laue_id"],
    )


# ---------------------------------------------------------------------------
# Misorientation functions (backward-compatible, now use the full lookup)
# ---------------------------------------------------------------------------

def misorientation_angle_sym(M: np.ndarray, symmetry: str = "triclinic") -> float:
    """
    Symmetry-reduced misorientation angle for a single (3, 3) relative rotation.

    Returns the disorientation angle (smallest equivalent rotation) in radians.
    """
    G = get_symmetry(symmetry)
    min_angle = np.inf
    for g1 in G:
        for g2 in G:
            Msym = g1 @ M @ g2.T
            trace = np.clip(np.trace(Msym), -1.0, 3.0)
            angle = np.arccos((trace - 1.0) / 2.0)
            if angle < min_angle:
                min_angle = angle
    return float(min_angle)


def misorientation_angles_batch(M_batch: np.ndarray,
                                symmetry: str = "triclinic") -> np.ndarray:
    """
    Vectorised symmetry-reduced misorientation angles.

    Parameters
    ----------
    M_batch  : (N, 3, 3) relative rotation matrices
    symmetry : crystal symmetry name

    Returns
    -------
    angles : (N,) disorientation angles in radians
    """
    G = get_symmetry(symmetry)
    N = len(M_batch)
    min_angles = np.full(N, np.inf)
    for g1 in G:
        for g2 in G:
            Msym = np.einsum("ij,njk,lk->nil", g1, M_batch, g2)
            trace = np.einsum("nii->n", Msym)
            angles = np.arccos(np.clip((trace - 1.0) / 2.0, -1.0, 1.0))
            min_angles = np.minimum(min_angles, angles)
    return min_angles


# ---------------------------------------------------------------------------
# CSL misorientations for cubic symmetry (Rodrigues-Frank space)
# Each entry: (angle_deg, axis_uvw)
# ---------------------------------------------------------------------------

CSL_CUBIC: dict[int, tuple[float, list[int]]] = {
    3:  (60.00, [1, 1, 1]),
    5:  (36.87, [1, 0, 0]),
    7:  (38.21, [1, 1, 1]),
    9:  (38.94, [1, 1, 0]),
    11: (50.48, [1, 1, 0]),
    13: (22.62, [1, 0, 0]),
    15: (48.19, [2, 1, 0]),
    17: (28.07, [1, 0, 0]),
    19: (26.53, [1, 1, 0]),
    21: (21.79, [1, 1, 1]),
    23: (40.45, [3, 1, 1]),
    25: (16.26, [1, 0, 0]),
    27: (31.59, [1, 1, 0]),
    29: (43.60, [1, 0, 0]),
    31: (17.90, [1, 1, 1]),
}


def csl_rotation_matrix(sigma: int) -> np.ndarray | None:
    """Return the reference rotation matrix for a given Σ value (cubic)."""
    if sigma not in CSL_CUBIC:
        return None
    angle_deg, axis = CSL_CUBIC[sigma]
    angle = np.radians(angle_deg)
    ax = np.array(axis, dtype=float)
    ax /= np.linalg.norm(ax)
    x, y, z = ax
    c, s = np.cos(angle), np.sin(angle)
    t = 1 - c
    R = np.array([
        [t*x*x + c,   t*x*y - s*z, t*x*z + s*y],
        [t*x*y + s*z, t*y*y + c,   t*y*z - s*x],
        [t*x*z - s*y, t*y*z + s*x, t*z*z + c  ],
    ])
    return R
