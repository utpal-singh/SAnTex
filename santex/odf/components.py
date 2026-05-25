"""
Standard texture components (named orientations) for common crystal systems.

Euler angles in degrees (Bunge ZXZ convention: φ1, Φ, φ2).
Each component may have several symmetry-equivalent variants.
"""

from __future__ import annotations
import numpy as np

# ---------------------------------------------------------------------------
# Cubic (m-3m) — rolling/deformation textures
# ---------------------------------------------------------------------------

CUBIC_COMPONENTS: dict[str, dict] = {
    "Cube":          {"euler": (0,   0,   0),  "description": "{001}<100>"},
    "RotatedCube":   {"euler": (45,  0,   0),  "description": "{001}<110>"},
    "Goss":          {"euler": (0,   45,  0),  "description": "{110}<001>"},
    "RotatedGoss":   {"euler": (90,  45,  0),  "description": "{110}<110>"},
    "Brass":         {"euler": (35,  45,  0),  "description": "{110}<112>"},
    "Copper":        {"euler": (90,  35,  45), "description": "{112}<111>"},
    "S":             {"euler": (59,  37,  63), "description": "{123}<634>"},
    "Dillamore":     {"euler": (90,  27,  45), "description": "{4 4 11}<11 11 8>"},
    "Taylor":        {"euler": (90,  27,  45), "description": "{4 4 11}<11 11 8>"},
    "P":             {"euler": (70,  45,  0),  "description": "{011}<111>"},
    "Q":             {"euler": (58,  18,  0),  "description": "{013}<100>"},
    "R":             {"euler": (57,  29,  63), "description": "{124}<211>"},
    "F":             {"euler": (54,  65,  0),  "description": "{111}<110>"},
    "E":             {"euler": (0,   90,  45), "description": "{011}<100>"},
}

# ---------------------------------------------------------------------------
# Hexagonal (6/mmm) — α-Ti, Mg, Zr etc.
# ---------------------------------------------------------------------------

HEX_COMPONENTS: dict[str, dict] = {
    "Basal":     {"euler": (0,   0,  0),  "description": "{0001}<10-10>"},
    "Prismatic": {"euler": (90,  90, 0),  "description": "{10-10}<0001>"},
    "Pyramidal": {"euler": (0,   45, 0),  "description": "{10-11}<11-20>"},
}

# ---------------------------------------------------------------------------
# Orthorhombic (mmm) — olivine / forsterite LPO fabrics
# ---------------------------------------------------------------------------

OLIVINE_FABRICS: dict[str, dict] = {
    # Euler angles approximate; exact values depend on the reference frame
    "A-type":  {"euler": (0,  0,  0),  "description": "[100](010) — [100]||X, [010]||Y"},
    "B-type":  {"euler": (0,  90, 0),  "description": "[100](001) — [001]||X"},
    "C-type":  {"euler": (90, 0,  0),  "description": "[001](010) — [001]||X"},
    "D-type":  {"euler": (45, 0,  0),  "description": "[100](010) — diffuse"},
    "E-type":  {"euler": (0,  0,  0),  "description": "[100](001) — [001]||Y"},
}

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def get_components(crystal_symmetry: str) -> dict[str, dict]:
    """Return the component dictionary for a given crystal symmetry."""
    s = crystal_symmetry.lower()
    if "cubic" in s or "m-3m" in s or "m3m" in s:
        return CUBIC_COMPONENTS
    if "hex" in s or "6/mmm" in s:
        return HEX_COMPONENTS
    if "ortho" in s or "mmm" in s or "olivine" in s or "forsterite" in s:
        return OLIVINE_FABRICS
    return CUBIC_COMPONENTS   # default


def volume_fraction(euler_data_deg: np.ndarray,
                    component_euler_deg: tuple[float, float, float],
                    max_angle_deg: float = 15.0,
                    weights: np.ndarray | None = None) -> float:
    """
    Compute the volume fraction of orientations within max_angle_deg of a
    texture component (no crystal symmetry reduction — use the raw angle).

    Parameters
    ----------
    euler_data_deg : (N, 3) Euler angles in degrees
    component_euler_deg : (phi1, PHI, phi2) of the component in degrees
    max_angle_deg : inclusion cone half-angle (degrees)
    weights : per-orientation weight array, shape (N,); None = uniform

    Returns
    -------
    fraction : float in [0, 1]
    """
    from santex.odf.odf import _euler_to_quat
    import numpy as _np

    euler_rad = _np.radians(euler_data_deg)
    q_data  = _euler_to_quat(euler_rad[:, 0], euler_rad[:, 1], euler_rad[:, 2])  # (N, 4)
    comp_rad = _np.radians(component_euler_deg)
    q_comp  = _euler_to_quat(*comp_rad)  # (4,)

    # Misorientation angle without symmetry reduction (fast approximation)
    dot = np.abs(q_data @ q_comp)          # (N,)
    dot = np.clip(dot, 0.0, 1.0)
    omega_deg = np.degrees(2.0 * np.arccos(dot))  # (N,)

    within = omega_deg <= max_angle_deg
    if weights is None:
        return float(within.mean())
    return float(weights[within].sum() / weights.sum())


def component_table(euler_data_deg: np.ndarray,
                    crystal_symmetry: str = "cubic",
                    max_angle_deg: float = 15.0,
                    weights: np.ndarray | None = None) -> list[dict]:
    """
    Compute volume fractions for all standard components.

    Returns
    -------
    rows : list of dicts with keys 'name', 'euler', 'description', 'fraction'
    """
    comps = get_components(crystal_symmetry)
    rows = []
    for name, info in comps.items():
        frac = volume_fraction(
            euler_data_deg, info["euler"], max_angle_deg, weights
        )
        rows.append({
            "name":        name,
            "euler":       info["euler"],
            "description": info["description"],
            "fraction":    frac,
        })
    return rows
