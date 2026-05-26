"""
Internal IPF colouring — pure NumPy, no orix dependency.

Architecture
------------
1.  Alias map — every HM / Schoenflies / short name → one of 11
    canonical Laue-group keys (C1, Ci, C2h, D2h, C4h, D4h, S6, D3d,
    C6h, D6h, Th, Oh).
2.  Space-group → Laue-group converter (`spacegroup_to_laue`).
    Accepts ITA space-group numbers 1–230 from CTF "Symmetry" fields.
3.  Elementary rotation-matrix generators; group-closure algorithm;
    per-group proper-rotation subgroup cache.
4.  Fundamental-sector data: corner poles, TSL corner colours, HSV
    sector phi range, and an internal reference direction for folding.
5.  Batch direction-folding via einsum + argmax.
6.  Two colour keys (matching MTEX):
      ``'tsl'``  – softmax corner interpolation (blue [001], green [0k1],
                   red [111]) — default; matches MTEX ipfTSLKey / ipfHKLKey.
      ``'hsv'``  – MTEX ipfHSVKey: white at pole, rainbow at boundary;
                   hue = azimuth within sector, saturation = sin(θ).
7.  `render_rgb_map` — scatter EBSD data → regular-grid RGBA uint8 image
    for `plotly.graph_objects.Image`.
8.  `make_colorkey_image` — renders the IPF colour-key triangle.

References
----------
* MTEX source: https://github.com/mtex-toolbox/mtex
* TSL color-key convention: TSL OIM Analysis User Manual
* Softmax IPF coloring: Nolze & Hielscher (2016) J. Appl. Cryst.
"""

from __future__ import annotations
import numpy as np

# ─────────────────────────────────────────────────────────────────────────────
# 1.  Space-group → Laue-group converter
# ─────────────────────────────────────────────────────────────────────────────

def spacegroup_to_laue(sg) -> str:
    """Return the canonical Laue-group key for an ITA space-group number.

    CTF files store the symmetry as an integer space-group number (1–230).
    Some older Channel-5 files store the Laue group number (1–11); both
    cases are handled here.

    Parameters
    ----------
    sg : int or str
        ITA space-group number (1–230), or Channel-5 Laue-group number
        (1–11), or a string representation of either.

    Returns
    -------
    str
        One of the 11 canonical Laue-group keys.
    """
    try:
        n = int(float(str(sg).split()[0]))
    except (ValueError, TypeError):
        return "D2h"

    # Channel-5 Laue group codes 1–11 (often stored directly):
    _LAUE_CODES = {
        1: "Ci", 2: "C2h", 3: "D2h", 4: "C4h", 5: "D4h",
        6: "S6",  7: "D3d", 8: "C6h", 9: "D6h", 10: "Th", 11: "Oh",
    }
    if 1 <= n <= 11:
        return _LAUE_CODES[n]

    # ITA space-group numbers 1–230:
    if   n == 1:            return "C1"
    elif n == 2:            return "Ci"
    elif 3  <= n <= 15:     return "C2h"
    elif 16 <= n <= 74:     return "D2h"
    elif 75 <= n <= 88:     return "C4h"
    elif 89 <= n <= 142:    return "D4h"
    elif 143 <= n <= 148:   return "S6"
    elif 149 <= n <= 167:   return "D3d"
    elif 168 <= n <= 176:   return "C6h"
    elif 177 <= n <= 194:   return "D6h"
    elif 195 <= n <= 206:   return "Th"
    elif 207 <= n <= 230:   return "Oh"
    else:                   return "D2h"


# ─────────────────────────────────────────────────────────────────────────────
# 2.  Canonical Laue-group keys and alias map
# ─────────────────────────────────────────────────────────────────────────────

LAUE_GROUPS: tuple[str, ...] = (
    "C1", "Ci",
    "C2h",
    "D2h",
    "C4h", "D4h",
    "S6",  "D3d",
    "C6h", "D6h",
    "Th",  "Oh",
)

_ALIAS: dict[str, str] = {
    # triclinic
    "1":     "C1",  "C1":  "C1",
    "-1":    "Ci",  "Ci":  "Ci",
    # monoclinic
    "2":     "C2h", "m":   "C2h", "2/m":  "C2h",
    "C2":    "C2h", "Cs":  "C2h", "C2h":  "C2h",
    # orthorhombic
    "222":   "D2h", "mm2": "D2h", "mmm":  "D2h",
    "D2":    "D2h", "C2v": "D2h", "D2h":  "D2h",
    # tetragonal
    "4":     "C4h", "-4":  "C4h", "4/m":  "C4h",
    "C4":    "C4h", "S4":  "C4h", "C4h":  "C4h",
    "422":   "D4h", "4mm": "D4h", "-42m": "D4h", "-4m2": "D4h",
    "4/mmm": "D4h",
    "D4":    "D4h", "C4v": "D4h", "D2d":  "D4h", "D4h":  "D4h",
    # trigonal
    "3":     "S6",  "-3":  "S6",
    "C3":    "S6",  "C3i": "S6",  "S6":   "S6",
    "32":    "D3d", "3m":  "D3d", "-3m":  "D3d",
    "321":   "D3d", "3m1": "D3d", "-3m1": "D3d",
    "312":   "D3d", "31m": "D3d", "-31m": "D3d",
    "D3":    "D3d", "C3v": "D3d", "D3d":  "D3d",
    # hexagonal
    "6":     "C6h", "-6":  "C6h", "6/m":  "C6h",
    "C6":    "C6h", "C3h": "C6h", "C6h":  "C6h",
    "622":   "D6h", "6mm": "D6h", "-62m": "D6h", "-6m2": "D6h",
    "6/mmm": "D6h",
    "D6":    "D6h", "C6v": "D6h", "D3h":  "D6h", "D6h":  "D6h",
    # cubic
    "23":    "Th",  "m-3": "Th",
    "T":     "Th",  "Th":  "Th",
    "432":   "Oh",  "-43m":"Oh",  "m-3m": "Oh",
    "O":     "Oh",  "Td":  "Oh",  "Oh":   "Oh",
}


def laue_group(name: str) -> str:
    """Return the canonical Laue-group key for *name*.

    Accepts Schoenflies (``'Oh'``, ``'D2h'``, …), Hermann–Mauguin short
    (``'m-3m'``, ``'mmm'``, …), and common variants (``'C3i'``, …).
    Also accepts space-group numbers passed as strings (``'62'``) — these
    are routed through :func:`spacegroup_to_laue`.
    Falls back to ``'D2h'`` for unknown input.
    """
    key = _ALIAS.get(name, _ALIAS.get(name.strip(), None))
    if key is not None:
        return key
    # Try as a space/Laue-group number
    try:
        return spacegroup_to_laue(int(float(name.strip())))
    except (ValueError, TypeError):
        return "D2h"


# ─────────────────────────────────────────────────────────────────────────────
# 3.  Elementary rotation-matrix generators
# ─────────────────────────────────────────────────────────────────────────────

_E = np.eye(3, dtype=np.float64)


def _Rz(deg: float) -> np.ndarray:
    a = np.deg2rad(deg)
    c, s = np.cos(a), np.sin(a)
    return np.array([[c, -s, 0.], [s, c, 0.], [0., 0., 1.]], dtype=np.float64)


def _C2n(n: np.ndarray) -> np.ndarray:
    nv = np.asarray(n, dtype=np.float64)
    nv = nv / np.linalg.norm(nv)
    return 2.0 * np.outer(nv, nv) - _E


_C2z  = np.diag([-1., -1.,  1.]).astype(np.float64)
_C2y  = np.diag([-1.,  1., -1.]).astype(np.float64)
_C2x  = np.diag([ 1., -1., -1.]).astype(np.float64)
_C4z  = np.array([[0., -1., 0.], [1., 0., 0.], [0., 0., 1.]], dtype=np.float64)
_C3z  = _Rz(120.)
_C6z  = _Rz( 60.)
# 3-fold about [111]: cycles x→y, y→z, z→x
_C3_111 = np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]], dtype=np.float64)


# ─────────────────────────────────────────────────────────────────────────────
# 4.  Group-closure algorithm + cache
# ─────────────────────────────────────────────────────────────────────────────

def _group_closure(
    generators: list[np.ndarray],
    tol: float = 1e-9,
    max_ops: int = 96,
) -> np.ndarray:
    ops: list[np.ndarray] = [_E.copy()]
    for g in generators:
        g_arr = np.asarray(g, dtype=np.float64)
        if not any(np.allclose(g_arr, o, atol=tol) for o in ops):
            ops.append(g_arr)
    prev_len = 0
    while len(ops) > prev_len:
        prev_len = len(ops)
        snapshot = list(ops)
        for a in snapshot:
            for b in snapshot:
                p = a @ b
                if not any(np.allclose(p, o, atol=tol) for o in ops):
                    ops.append(p)
                    if len(ops) >= max_ops:
                        return np.stack(ops, axis=0)
    return np.stack(ops, axis=0)


_SYM_CACHE: dict[str, np.ndarray] = {}


def get_sym_ops(name: str) -> np.ndarray:
    """Return (N, 3, 3) proper-rotation operators for the Laue group *name*."""
    key = laue_group(name)
    if key in _SYM_CACHE:
        return _SYM_CACHE[key]
    if   key == "C1":  ops = _E[np.newaxis].copy()
    elif key == "Ci":  ops = _E[np.newaxis].copy()
    elif key == "C2h": ops = _group_closure([_C2z])
    elif key == "D2h": ops = _group_closure([_C2z, _C2x])
    elif key == "C4h": ops = _group_closure([_C4z])
    elif key == "D4h": ops = _group_closure([_C4z, _C2x])
    elif key == "S6":  ops = _group_closure([_C3z])
    elif key == "D3d": ops = _group_closure([_C3z, _C2x])
    elif key == "C6h": ops = _group_closure([_C6z])
    elif key == "D6h": ops = _group_closure([_C6z, _C2x])
    elif key == "Th":  ops = _group_closure([_C2z, _C2x, _C3_111])
    elif key == "Oh":  ops = _group_closure([_C4z, _C3_111])
    else:              ops = _E[np.newaxis].copy()
    _SYM_CACHE[key] = ops
    return ops


# ─────────────────────────────────────────────────────────────────────────────
# 5.  Fundamental-sector data
# ─────────────────────────────────────────────────────────────────────────────

_s2 = np.sqrt(2.)
_s3 = np.sqrt(3.)

# TSL corner colours ─ [001]=blue, [0k1]/[k00]=green or red,  [111]=red
# HSV phi range ─ (phi_min, phi_max) in radians
#   phi = atan2(hy, hx) within the sector
#   H_norm = (phi - phi_min) / (phi_max - phi_min)  ∈ [0, 1]

_SECTOR: dict[str, dict] = {

    # ── triclinic ────────────────────────────────────────────────────────────
    "C1": {
        "corners": np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]]),
        "colors":  np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]]),
        "ref":     np.array([1., 1., 1.]) / _s3,
        "phi":     (0., 2*np.pi),
    },
    "Ci": {
        "corners": np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]]),
        "colors":  np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]]),
        "ref":     np.array([1., 1., 1.]) / _s3,
        "phi":     (0., 2*np.pi),
    },

    # ── monoclinic ───────────────────────────────────────────────────────────
    "C2h": {
        "corners": np.array([[0., 0., 1.], [1., 0., 0.]]),
        "colors":  np.array([[0., 0., 1.], [1., 0., 0.]]),
        "ref":     np.array([1., 0., 1.]) / _s2,
        "phi":     (0., np.pi),
    },

    # ── orthorhombic ─────────────────────────────────────────────────────────
    # Sector: z≥y≥x≥0 with folding reference (1,1,1)/√3.
    # HKL convention: [001]=blue, [010]=green, [100]=red.
    # phi range: [100] has phi=0°, [010] has phi=90°  →  0 to π/2.
    "D2h": {
        "corners": np.array([[0., 0., 1.], [0., 1., 0.], [1., 0., 0.]]),
        "colors":  np.array([[0., 0., 1.], [0., 1., 0.], [1., 0., 0.]]),
        "ref":     np.array([1., 1., 1.]) / _s3,
        "phi":     (0., np.pi / 2),
    },

    # ── tetragonal ───────────────────────────────────────────────────────────
    "C4h": {
        "corners": np.array([[0., 0., 1.], [1., 0., 0.]]),
        "colors":  np.array([[0., 0., 1.], [1., 0., 0.]]),
        "ref":     np.array([1., 0., 1.]) / _s2,
        "phi":     (0., np.pi / 2),
    },
    # Sector corners: [001], [100], [110]/√2
    # phi range: [100] (phi=0) to [110]/√2 (phi=45°)  →  0 to π/4.
    "D4h": {
        "corners": np.array([
            [0., 0., 1.],
            [1., 0., 0.],
            [1./_s2, 1./_s2, 0.],
        ]),
        "colors": np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]]),
        "ref":    np.array([1., 0.5, 2.]) / np.linalg.norm([1., 0.5, 2.]),
        "phi":    (0., np.pi / 4),
    },

    # ── trigonal ─────────────────────────────────────────────────────────────
    "S6": {
        "corners": np.array([[0., 0., 1.], [1., 0., 0.]]),
        "colors":  np.array([[0., 0., 1.], [1., 0., 0.]]),
        "ref":     np.array([1., 0., 1.]) / _s2,
        "phi":     (0., 2*np.pi / 3),
    },
    # Sector corners: [001], [100], [1/2, √3/2, 0]  (azimuth 0° → 60°)
    "D3d": {
        "corners": np.array([
            [0.,   0.,        1.],
            [1.,   0.,        0.],
            [0.5,  _s3 / 2.,  0.],
        ]),
        "colors": np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]]),
        "ref":    np.array([np.cos(np.deg2rad(30.)),
                             np.sin(np.deg2rad(30.)), 1.]),
        "phi":    (0., np.pi / 3),
    },

    # ── hexagonal ────────────────────────────────────────────────────────────
    "C6h": {
        "corners": np.array([[0., 0., 1.], [1., 0., 0.]]),
        "colors":  np.array([[0., 0., 1.], [1., 0., 0.]]),
        "ref":     np.array([1., 0., 1.]) / _s2,
        "phi":     (0., np.pi / 3),
    },
    # Sector corners: [001], [100], [√3/2, 1/2, 0]  (azimuth 0° → 30°)
    "D6h": {
        "corners": np.array([
            [0.,       0.,   1.],
            [1.,       0.,   0.],
            [_s3/2.,   0.5,  0.],
        ]),
        "colors": np.array([[0., 0., 1.], [1., 0., 0.], [0., 1., 0.]]),
        "ref":    np.array([np.cos(np.deg2rad(15.)),
                             np.sin(np.deg2rad(15.)), 1.]),
        "phi":    (0., np.pi / 6),
    },

    # ── cubic ────────────────────────────────────────────────────────────────
    "Th": {
        "corners": np.array([
            [0.,      0.,      1.],
            [0.,      1./_s2,  1./_s2],
            [1./_s3,  1./_s3,  1./_s3],
        ]),
        "colors": np.array([[0., 0., 1.], [0., 1., 0.], [1., 0., 0.]]),
        "ref":    np.array([1., 2., 3.]) / np.linalg.norm([1., 2., 3.]),
        "phi":    (np.pi / 4, np.pi / 2),
    },
    # Cubic m-3m (HKL convention): [001]=blue, [011]/√2=green, [111]/√3=red
    # Sector z≥y≥x≥0; phi: [111] (phi=45°) → [011] (phi=90°)  →  π/4 to π/2.
    "Oh": {
        "corners": np.array([
            [0.,      0.,      1.],
            [0.,      1./_s2,  1./_s2],
            [1./_s3,  1./_s3,  1./_s3],
        ]),
        "colors": np.array([[0., 0., 1.], [0., 1., 0.], [1., 0., 0.]]),
        "ref":    np.array([1., 2., 3.]) / np.linalg.norm([1., 2., 3.]),
        "phi":    (np.pi / 4, np.pi / 2),
    },
}

# Normalise corners and reference in-place
for _key, _sec in _SECTOR.items():
    _c = _sec["corners"].astype(np.float64)
    _n = np.linalg.norm(_c, axis=1, keepdims=True)
    _sec["corners"] = np.where(_n > 0, _c / _n, _c)
    _r = np.asarray(_sec["ref"], dtype=np.float64)
    _sec["ref"] = _r / np.linalg.norm(_r)
del _key, _sec, _c, _n, _r


# ─────────────────────────────────────────────────────────────────────────────
# Sector corner labels  (used for IPF colour-key annotation)
# ─────────────────────────────────────────────────────────────────────────────
# Each entry is a list of (xyz_direction, label_text) pairs in order:
#   first entry = pole direction ([001] / [0001]) — displayed near the apex
#   remaining  = equatorial vertices

_SECTOR_CORNER_LABELS: dict[str, list[tuple[list[float], str]]] = {
    "C1":  [([0, 0, 1], "[001]")],
    "Ci":  [([0, 0, 1], "[001]"), ([1, 0, 0], "[100]")],
    "C2h": [([0, 0, 1], "[001]"), ([1, 0, 0], "[100]")],
    "D2h": [([0, 0, 1], "[001]"), ([0, 1, 0], "[010]"), ([1, 0, 0], "[100]")],
    "C4h": [([0, 0, 1], "[001]"), ([1, 0, 0], "[100]")],
    "D4h": [([0, 0, 1], "[001]"), ([1, 0, 0], "[100]"),
            ([1/_s2, 1/_s2, 0], "[110]")],
    "S6":  [([0, 0, 1], "[0001]"), ([1, 0, 0], "[10̐10]")],
    "D3d": [([0, 0, 1], "[0001]"), ([1, 0, 0], "[10̐10]"),
            ([0.5, _s3/2., 0], "[11̐2̄0]")],
    "C6h": [([0, 0, 1], "[0001]"), ([1, 0, 0], "[10̐10]")],
    "D6h": [([0, 0, 1], "[0001]"), ([1, 0, 0], "[10̐10]"),
            ([_s3/2., 0.5, 0], "[11̐2̄0]")],
    "Th":  [([0, 0, 1], "[001]"), ([0, 1/_s2, 1/_s2], "[011]"),
            ([1/_s3, 1/_s3, 1/_s3], "[111]")],
    "Oh":  [([0, 0, 1], "[001]"), ([0, 1/_s2, 1/_s2], "[011]"),
            ([1/_s3, 1/_s3, 1/_s3], "[111]")],
}


def sector_corner_pixels(
    sector_key: str,
    n: int = 256,
) -> list[tuple[str, int, int, str, str]]:
    """Return annotation data for each fundamental-sector corner.

    Pixel coordinates correspond to the output of
    ``make_colorkey_image(sector_key, n=n)``.

    Parameters
    ----------
    sector_key : canonical Laue-group key (e.g. ``'Oh'``, ``'D2h'``)
    n          : image side length (pixels) — must match ``make_colorkey_image``

    Returns
    -------
    list of (label, col_px, row_px, xanchor, yanchor)
        Plotly-compatible anchor strings to position the label *away* from
        the centre of the triangle.
    """
    corners_info = _SECTOR_CORNER_LABELS.get(sector_key, [([0, 0, 1], "[001]")])
    cx = n // 2   # image centre column
    cy = n // 2   # image centre row

    result: list[tuple[str, int, int, str, str]] = []
    for xyz, label in corners_info:
        v   = np.asarray(xyz, dtype=np.float64)
        nrm = np.linalg.norm(v)
        if nrm < 1e-10:                      # degenerate
            result.append((label, cx, cy, "center", "bottom"))
            continue

        v   = v / nrm
        z   = float(np.clip(v[2], -1.0, 1.0))
        # Lambert equal-area: r = √2 · sin(θ/2)
        theta  = np.arccos(z)
        r      = np.sqrt(2.0) * np.sin(theta / 2.0)
        phi    = np.arctan2(v[1], v[0])
        xi     = r * np.cos(phi)
        yi     = r * np.sin(phi)

        # Convert Lambert grid (−1 … 1) → pixel index (0 … n−1)
        col_px = int(round((xi + 1.0) / 2.0 * (n - 1)))
        row_px = int(round((yi + 1.0) / 2.0 * (n - 1)))

        # Anchor: push label AWAY from the image centre
        dx = col_px - cx
        dy = row_px - cy
        thresh = n // 6

        if abs(dx) <= thresh and abs(dy) <= thresh:
            # Near centre — always put text above
            xanchor, yanchor = "center", "bottom"
        else:
            xanchor = "left"  if dx >  thresh else (
                      "right" if dx < -thresh else "center")
            yanchor = "top"   if dy >  thresh else (
                      "bottom" if dy < -thresh else "middle")

        result.append((label, col_px, row_px, xanchor, yanchor))

    return result


# ─────────────────────────────────────────────────────────────────────────────
# 6a.  Batch direction folding
# ─────────────────────────────────────────────────────────────────────────────

def fold_directions(V: np.ndarray, sym_key: str) -> np.ndarray:
    """Map each direction in *V* to the canonical fundamental-sector representative.

    Parameters
    ----------
    V       : (N, 3) unit vectors (crystal frame)
    sym_key : canonical Laue-group key from :func:`laue_group`

    Returns
    -------
    (N, 3) unit vectors in the fundamental sector.
    """
    ops = get_sym_ops(sym_key)
    ref = _SECTOR.get(sym_key, _SECTOR["D2h"])["ref"]

    equiv  = np.einsum("jkl,il->ijk", ops, V)          # (N, M, 3)
    equiv  = np.concatenate([equiv, -equiv], axis=1)    # (N, 2M, 3)
    scores = np.einsum("ijk,k->ij", equiv, ref)          # (N, 2M)
    best   = np.argmax(scores, axis=1)
    N      = len(V)
    canon  = equiv[np.arange(N), best]
    norms  = np.linalg.norm(canon, axis=1, keepdims=True)
    return np.where(norms > 1e-15, canon / norms, canon)


# ─────────────────────────────────────────────────────────────────────────────
# 6b.  TSL colour key  (softmax corner interpolation)
# ─────────────────────────────────────────────────────────────────────────────

def _ipf_tsl_color(v_canonical: np.ndarray, sector_key: str,
                   T: float = 50.0) -> np.ndarray:
    """TSL / HKL softmax colour key.

    Corners [001], [0k1], [111] mapped to blue, green, red; intermediate
    directions interpolated via softmax on arccos distance.

    Parameters
    ----------
    v_canonical : (N, 3) unit directions in fundamental sector
    sector_key  : canonical Laue-group key
    T           : softmax sharpness (50 ≈ MTEX default)

    Returns
    -------
    (N, 3) RGB in [0, 1]
    """
    sec     = _SECTOR.get(sector_key, _SECTOR["D2h"])
    corners = sec["corners"]
    colors  = sec["colors"]

    cos_ang = np.clip(v_canonical @ corners.T, -1.0, 1.0)
    angles  = np.arccos(cos_ang)
    weights = np.exp(-T * angles ** 2)
    weights /= weights.sum(axis=1, keepdims=True) + 1e-30
    return np.clip(weights @ colors, 0.0, 1.0)


# ─────────────────────────────────────────────────────────────────────────────
# 6c.  HSV colour key  (MTEX ipfHSVKey — white at pole, rainbow at boundary)
# ─────────────────────────────────────────────────────────────────────────────

def _hsv_to_rgb_batch(H: np.ndarray, S: np.ndarray,
                      V: np.ndarray) -> np.ndarray:
    """Vectorised HSV → RGB.  All inputs/outputs are (N,) arrays in [0, 1]."""
    H6  = H * 6.0
    I   = np.floor(H6).astype(int) % 6
    F   = H6 - np.floor(H6)
    P   = V * (1.0 - S)
    Q   = V * (1.0 - S * F)
    Tv  = V * (1.0 - S * (1.0 - F))
    N   = len(H)
    rgb = np.zeros((N, 3), dtype=np.float64)
    for k, (r, g, b) in enumerate([
        (V, Tv, P), (Q, V,  P),
        (P, V,  Tv),(P, Q,  V),
        (Tv,P,  V), (V, P,  Q),
    ]):
        m = (I == k)
        rgb[m, 0] = r[m]; rgb[m, 1] = g[m]; rgb[m, 2] = b[m]
    return rgb


def _ipf_hsv_color(v_canonical: np.ndarray, sector_key: str,
                   hue_scale: float = 2.0 / 3.0) -> np.ndarray:
    """MTEX ``ipfHSVKey`` colour mapping.

    * **Hue** (H) — azimuthal angle within the sector, normalised to
      [0, *hue_scale*].  The default *hue_scale = 2/3* maps the sector
      across the red→green→blue portion of the HSV wheel, matching MTEX.
    * **Saturation** (S) — sin(θ) where θ is the polar angle from the
      main pole ([001]); S = 0 at the pole (white), S = 1 at the equator.
    * **Value** = 1 (full brightness).

    Parameters
    ----------
    v_canonical : (N, 3) unit directions in fundamental sector
    sector_key  : canonical Laue-group key
    hue_scale   : fraction of the HSV hue wheel used (default 2/3)

    Returns
    -------
    (N, 3) RGB in [0, 1]
    """
    sec                 = _SECTOR.get(sector_key, _SECTOR["D2h"])
    phi_min, phi_max    = sec["phi"]
    phi_span            = phi_max - phi_min

    hx = v_canonical[:, 0]
    hy = v_canonical[:, 1]
    hz = np.clip(v_canonical[:, 2], -1.0, 1.0)

    # Azimuthal angle [0, 2π)
    phi = np.arctan2(hy, hx) % (2.0 * np.pi)

    # Normalised hue within sector [0, hue_scale]
    if phi_span > 1e-10:
        H = np.clip((phi - phi_min) / phi_span, 0.0, 1.0) * hue_scale
    else:
        H = np.zeros(len(v_canonical))

    # Saturation = sin(polar angle from z-pole)
    S = np.clip(np.sqrt(hx**2 + hy**2), 0.0, 1.0)

    V = np.ones(len(v_canonical))
    return _hsv_to_rgb_batch(H, S, V)


# ─────────────────────────────────────────────────────────────────────────────
# 7.  Private rotation utilities
# ─────────────────────────────────────────────────────────────────────────────

def _euler_to_rotation_batch(
    phi1: np.ndarray, PHI: np.ndarray, phi2: np.ndarray
) -> np.ndarray:
    """Vectorised Bunge ZXZ → rotation matrix (N, 3, 3)."""
    c1, s1 = np.cos(phi1), np.sin(phi1)
    c,  s  = np.cos(PHI),  np.sin(PHI)
    c2, s2 = np.cos(phi2), np.sin(phi2)
    N = len(phi1)
    R = np.zeros((N, 3, 3), dtype=np.float64)
    R[:, 0, 0] =  c1*c2 - s1*s2*c;  R[:, 0, 1] = -c1*s2 - s1*c2*c
    R[:, 0, 2] =  s1*s
    R[:, 1, 0] =  s1*c2 + c1*s2*c;  R[:, 1, 1] = -s1*s2 + c1*c2*c
    R[:, 1, 2] = -c1*s
    R[:, 2, 0] =  s2*s;              R[:, 2, 1] =  c2*s
    R[:, 2, 2] =  c
    return R


# ─────────────────────────────────────────────────────────────────────────────
# 8.  Public API — ipf_map_colors
# ─────────────────────────────────────────────────────────────────────────────

def ipf_map_colors(
    data,
    phase_index: int = 1,
    direction: str = "ND",
    crystal_symmetry: str = "D2h",
    color_key: str = "tsl",
) -> dict | None:
    """Compute per-pixel IPF colour for *phase_index* — no orix dependency.

    Parameters
    ----------
    data             : EBSD DataFrame (Phase, X, Y, Euler1, Euler2, Euler3).
    phase_index      : phase number (Phase == 0 pixels are excluded).
    direction        : sample reference direction —
                       ``'ND'`` (Z), ``'RD'`` (X), ``'TD'`` (Y).
    crystal_symmetry : symmetry name **or** a space-group / Laue-group number
                       as a string (``'62'``, ``'11'``, ``'Oh'``, ``'mmm'``…).
                       Pass ``'auto'`` to let the caller handle auto-detection
                       (the string is then resolved via :func:`laue_group`).
    color_key        : ``'tsl'`` (default) — blue/green/red corner colours;
                       ``'hsv'`` — MTEX ``ipfHSVKey``: white at pole,
                       rainbow at boundary.

    Returns
    -------
    dict with keys ``x``, ``y``, ``r``, ``g``, ``b`` (float 0–1),
    ``rgb_hex`` (list[str]), ``sym_key`` (str Laue key used).
    Returns ``None`` if no indexed pixels for *phase_index*.
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

    phi1 = np.deg2rad(phase_data["Euler1"].to_numpy(dtype=float))
    PHI  = np.deg2rad(phase_data["Euler2"].to_numpy(dtype=float))
    phi2 = np.deg2rad(phase_data["Euler3"].to_numpy(dtype=float))
    R = _euler_to_rotation_batch(phi1, PHI, phi2)            # (N, 3, 3)

    # Crystal direction that aligns with sample_dir:
    #   d_c = R^T @ d_s  (passive convention: R maps crystal → sample)
    # einsum "nji,j->ni" = R[n]^T @ d_s
    crystal_dirs = np.einsum("nji,j->ni", R, sample_dir)     # (N, 3)
    norms        = np.linalg.norm(crystal_dirs, axis=1, keepdims=True)
    crystal_dirs = np.where(norms > 1e-15, crystal_dirs / norms, crystal_dirs)

    sym_key   = laue_group(crystal_symmetry)
    canonical = fold_directions(crystal_dirs, sym_key)

    if color_key == "hsv":
        rgb = _ipf_hsv_color(canonical, sym_key)
    else:
        rgb = _ipf_tsl_color(canonical, sym_key)

    hex_colors = [
        "#{:02x}{:02x}{:02x}".format(int(r*255), int(g*255), int(b*255))
        for r, g, b in rgb
    ]
    return {
        "x":       phase_data["X"].to_numpy(),
        "y":       phase_data["Y"].to_numpy(),
        "r":       rgb[:, 0],
        "g":       rgb[:, 1],
        "b":       rgb[:, 2],
        "rgb_hex": hex_colors,
        "sym_key": sym_key,
    }


# ─────────────────────────────────────────────────────────────────────────────
# 9.  Pixel-grid renderer — scatter data → regular RGBA image for go.Image
# ─────────────────────────────────────────────────────────────────────────────

def _estimate_step(x: np.ndarray, y: np.ndarray) -> float:
    """Estimate EBSD step size from median NN distance along x."""
    ux = np.unique(np.round(np.sort(x), 6))
    if len(ux) < 2:
        uy = np.unique(np.round(np.sort(y), 6))
        if len(uy) < 2:
            return 1.0
        return float(np.median(np.diff(uy)))
    return float(np.median(np.diff(ux)))


def render_rgb_map(
    x: np.ndarray, y: np.ndarray,
    r: np.ndarray, g: np.ndarray, b: np.ndarray,
    bg_gray: float = 0.15,
) -> tuple[np.ndarray, float, float, float, float]:
    """Convert scattered IPF RGB values to a regular-grid RGBA uint8 image.

    Parameters
    ----------
    x, y      : pixel centre coordinates (µm)
    r, g, b   : float RGB values in [0, 1]
    bg_gray   : grey level for unoccupied pixels (default 0.15)

    Returns
    -------
    img       : (H, W, 4) uint8 RGBA  (origin = top-left, y-down)
    x0, x1   : x extent
    y0, y1   : y extent  (y0 < y1; plotly Image uses y0 as the top)
    """
    step = _estimate_step(x, y)
    if step <= 0:
        step = 1.0

    x_min, x_max = float(x.min()), float(x.max())
    y_min, y_max = float(y.min()), float(y.max())

    nx = max(1, int(round((x_max - x_min) / step)) + 1)
    ny = max(1, int(round((y_max - y_min) / step)) + 1)

    bg = int(bg_gray * 255)
    img = np.full((ny, nx, 4), bg, dtype=np.uint8)
    img[:, :, 3] = 255   # full alpha everywhere

    xi = np.clip(np.round((x - x_min) / step).astype(int), 0, nx - 1)
    # EBSD y grows downward in the scan; keep that convention (row 0 = y_min)
    yi = np.clip(np.round((y - y_min) / step).astype(int), 0, ny - 1)

    img[yi, xi, 0] = (np.clip(r, 0., 1.) * 255).astype(np.uint8)
    img[yi, xi, 1] = (np.clip(g, 0., 1.) * 255).astype(np.uint8)
    img[yi, xi, 2] = (np.clip(b, 0., 1.) * 255).astype(np.uint8)
    img[yi, xi, 3] = 255

    return img, x_min, x_max, y_min, y_max


# ─────────────────────────────────────────────────────────────────────────────
# 10.  IPF colour-key triangle image  (for legend display)
# ─────────────────────────────────────────────────────────────────────────────

def make_colorkey_image(
    sector_key: str = "Oh",
    color_key:  str = "tsl",
    n: int = 256,
    T: float = 50.0,
) -> np.ndarray:
    """Render the IPF colour-key triangle as an (n, n, 4) RGBA uint8 array.

    The image uses Lambert equal-area projection (upper hemisphere).  Only
    pixels inside the fundamental sector are coloured; outside pixels are
    transparent (alpha=0).

    Parameters
    ----------
    sector_key : canonical Laue-group key
    color_key  : ``'tsl'`` or ``'hsv'``
    n          : image width = height in pixels
    T          : softmax sharpness (for TSL key)

    Returns
    -------
    (n, n, 4) uint8 RGBA
    """
    sec = _SECTOR.get(sector_key, _SECTOR["D2h"])
    ref = sec["ref"]

    # Build a regular Cartesian grid on the Lambert disk
    lin   = np.linspace(-1.0, 1.0, n)
    XI, YI = np.meshgrid(lin, lin)
    r2    = XI**2 + YI**2
    valid = r2 <= 1.0

    r_safe  = np.sqrt(np.where(valid, r2, 0.0))
    # Lambert inverse: r = √2·sin(θ/2)  →  θ = 2·arcsin(r/√2)
    theta   = 2.0 * np.arcsin(np.clip(r_safe / np.sqrt(2.0), 0.0, 1.0))
    phi     = np.arctan2(YI, XI) % (2.0 * np.pi)

    hx = np.sin(theta) * np.cos(phi)
    hy = np.sin(theta) * np.sin(phi)
    hz = np.cos(theta)
    V  = np.stack([hx.ravel(), hy.ravel(), hz.ravel()], axis=1)  # (n², 3)

    # Fold to sector
    ops    = get_sym_ops(sector_key)
    equiv  = np.einsum("jkl,il->ijk", ops, V)
    equiv  = np.concatenate([equiv, -equiv], axis=1)
    scores = np.einsum("ijk,k->ij", equiv, ref)
    best   = np.argmax(scores, axis=1)
    N_flat = len(V)
    canon  = equiv[np.arange(N_flat), best]
    norms  = np.linalg.norm(canon, axis=1, keepdims=True)
    canon  = np.where(norms > 1e-15, canon / norms, canon)

    # Only show pixels inside the fundamental sector
    # (dot product with ref must be the highest among all equivalents)
    score_best = scores[np.arange(N_flat), best]
    score_2nd  = np.partition(scores, -2, axis=1)[:, -2]
    in_sector  = (score_best - score_2nd) > 1e-4   # strictly best

    if color_key == "hsv":
        rgb_flat = _ipf_hsv_color(canon, sector_key)
    else:
        rgb_flat = _ipf_tsl_color(canon, sector_key, T=T)

    img = np.zeros((n, n, 4), dtype=np.uint8)
    mask = valid.ravel() & in_sector

    rgb8     = (np.clip(rgb_flat, 0., 1.) * 255).astype(np.uint8)
    img_flat = img.reshape(-1, 4)
    img_flat[mask, 0] = rgb8[mask, 0]
    img_flat[mask, 1] = rgb8[mask, 1]
    img_flat[mask, 2] = rgb8[mask, 2]
    img_flat[mask, 3] = 255

    return img.reshape(n, n, 4)
