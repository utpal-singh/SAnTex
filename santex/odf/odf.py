"""
ODF (Orientation Distribution Function) — kernel density estimation on SO(3).

The ODF is computed by distributing de la Vallée-Poussin kernels centred on
each measured orientation over the orientation space, accounting for crystal
symmetry.  Results are expressed in multiples of a random distribution (m.r.d.).

Reference
---------
Hielscher & Schaeben (2008) "A novel pole figure inversion method: specification
of the MTEX algorithm", J. Appl. Cryst. 41, 1024–1037.
"""

from __future__ import annotations
import numpy as np
from functools import cached_property


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _euler_to_quat(phi1, PHI, phi2):
    """Vectorised Bunge ZXZ Euler (radians) → quaternion [w,x,y,z]."""
    c1, s1 = np.cos(phi1 / 2), np.sin(phi1 / 2)
    c2, s2 = np.cos(PHI  / 2), np.sin(PHI  / 2)
    c3, s3 = np.cos(phi2 / 2), np.sin(phi2 / 2)
    w = c1 * c2 * c3 - s1 * c2 * s3
    x = c1 * s2 * c3 + s1 * s2 * s3
    y = s1 * s2 * c3 - c1 * s2 * s3
    z = s1 * c2 * c3 + c1 * c2 * s3
    return np.stack([w, x, y, z], axis=-1)


def _quat_mul_batch(q1, q2):
    """
    q1: (..., 4), q2: (..., 4) → (..., 4) quaternion product.
    Broadcasts over leading dimensions.
    """
    w1, x1, y1, z1 = q1[..., 0], q1[..., 1], q1[..., 2], q1[..., 3]
    w2, x2, y2, z2 = q2[..., 0], q2[..., 1], q2[..., 2], q2[..., 3]
    return np.stack([
        w1*w2 - x1*x2 - y1*y2 - z1*z2,
        w1*x2 + x1*w2 + y1*z2 - z1*y2,
        w1*y2 - x1*z2 + y1*w2 + z1*x2,
        w1*z2 + x1*y2 - y1*x2 + z1*w2,
    ], axis=-1)


def _rot_to_quat(R):
    """Single 3×3 rotation matrix → quaternion [w,x,y,z]."""
    trace = R[0, 0] + R[1, 1] + R[2, 2]
    if trace > 0:
        s = 0.5 / np.sqrt(trace + 1.0)
        return np.array([0.25 / s,
                         (R[2, 1] - R[1, 2]) * s,
                         (R[0, 2] - R[2, 0]) * s,
                         (R[1, 0] - R[0, 1]) * s])
    elif R[0, 0] > R[1, 1] and R[0, 0] > R[2, 2]:
        s = 2.0 * np.sqrt(1.0 + R[0, 0] - R[1, 1] - R[2, 2])
        return np.array([(R[2, 1] - R[1, 2]) / s, 0.25 * s,
                         (R[0, 1] + R[1, 0]) / s, (R[0, 2] + R[2, 0]) / s])
    elif R[1, 1] > R[2, 2]:
        s = 2.0 * np.sqrt(1.0 + R[1, 1] - R[0, 0] - R[2, 2])
        return np.array([(R[0, 2] - R[2, 0]) / s, (R[0, 1] + R[1, 0]) / s,
                         0.25 * s, (R[1, 2] + R[2, 1]) / s])
    else:
        s = 2.0 * np.sqrt(1.0 + R[2, 2] - R[0, 0] - R[1, 1])
        return np.array([(R[1, 0] - R[0, 1]) / s, (R[0, 2] + R[2, 0]) / s,
                         (R[1, 2] + R[2, 1]) / s, 0.25 * s])


def _symmetry_quaternions(crystal_symmetry: str) -> np.ndarray:
    """Return (S, 4) quaternion array for all proper crystal symmetry ops."""
    from santex.grains.symmetry import SYMMETRY
    ops = SYMMETRY.get(crystal_symmetry.lower(), SYMMETRY["triclinic"])
    return np.array([_rot_to_quat(R) for R in ops], dtype=float)  # (S, 4)


# ---------------------------------------------------------------------------
# ODF Euler-angle grid
# ---------------------------------------------------------------------------

def _euler_grid(phi2_values_deg, phi1_range_deg=(0, 360),
                PHI_range_deg=(0, 90), resolution_deg=5.0):
    """
    Build a 3-column (phi1, PHI, phi2) Euler angle grid in degrees.

    phi2 is fixed per section; phi1 and PHI are swept at given resolution.
    """
    phi1 = np.arange(phi2_values_deg[0], phi1_range_deg[1] + 1e-9, resolution_deg)
    PHI  = np.arange(PHI_range_deg[0],    PHI_range_deg[1]  + 1e-9, resolution_deg)
    grids = {}
    for phi2 in phi2_values_deg:
        p1, P = np.meshgrid(phi1, PHI)
        p2 = np.full_like(p1, phi2)
        grid = np.stack([p1.ravel(), P.ravel(), p2.ravel()], axis=1)
        grids[phi2] = (p1, P, grid)
    return grids


# ---------------------------------------------------------------------------
# ODF class
# ---------------------------------------------------------------------------

class ODF:
    """
    Orientation Distribution Function estimated by kernel density estimation.

    Parameters
    ----------
    euler_deg        : (N, 3) Euler angles in degrees (Bunge ZXZ)
    weights          : (N,) optional per-orientation weights
    crystal_symmetry : crystal symmetry identifier string
    sample_symmetry  : 'none' (triclinic), 'orthorhombic', 'monoclinic'
    halfwidth_deg    : FWHM half-width of the de la Vallée-Poussin kernel
    max_orientations : subsample to this many orientations for speed
    """

    def __init__(self,
                 euler_deg: np.ndarray,
                 weights: np.ndarray | None = None,
                 crystal_symmetry: str = "orthorhombic",
                 sample_symmetry: str = "none",
                 halfwidth_deg: float = 10.0,
                 max_orientations: int = 5000):

        self.crystal_symmetry = crystal_symmetry
        self.sample_symmetry  = sample_symmetry
        self.halfwidth_deg    = halfwidth_deg

        euler_deg = np.asarray(euler_deg, dtype=float)
        if weights is not None:
            weights = np.asarray(weights, dtype=float)
            weights = weights / weights.sum()

        # Sub-sample if too many orientations
        N = len(euler_deg)
        if N > max_orientations:
            rng = np.random.default_rng(42)
            idx = rng.choice(N, max_orientations, replace=False)
            euler_deg = euler_deg[idx]
            weights = weights[idx] / weights[idx].sum() if weights is not None else None
            N = max_orientations

        self._euler_deg = euler_deg
        self._weights   = weights
        self._N         = N
        self._kappa     = _halfwidth_to_kappa(halfwidth_deg)
        self._precompute()

    # ------------------------------------------------------------------

    def _precompute(self):
        """Pre-compute data quaternions × all symmetry equivalents."""
        from .kernel import halfwidth_to_kappa
        self._kappa = halfwidth_to_kappa(self.halfwidth_deg)

        euler_rad = np.radians(self._euler_deg)
        q_data = _euler_to_quat(euler_rad[:, 0], euler_rad[:, 1], euler_rad[:, 2])
        # Ensure positive w hemisphere (canonical form)
        q_data = q_data * np.sign(q_data[:, 0:1] + 1e-30)

        q_sym = _symmetry_quaternions(self.crystal_symmetry)  # (S, 4)
        S = len(q_sym)

        # q_sym_data[n, s] = q_data[n] * q_sym[s]
        # Shape: (N, S, 4)
        q_ds = _quat_mul_batch(
            q_data[:, None, :],   # (N, 1, 4)
            q_sym[None, :, :]     # (1, S, 4)
        )   # (N, S, 4)
        # Flatten to (N*S, 4) for fast matrix multiply
        self._q_sym_data = q_ds.reshape(-1, 4)   # (N*S, 4)
        self._S = S

    # ------------------------------------------------------------------
    # Core evaluation
    # ------------------------------------------------------------------

    def evaluate(self, euler_grid_deg: np.ndarray,
                 batch_size: int = 200) -> np.ndarray:
        """
        Evaluate ODF (in m.r.d.) at a set of Euler angle grid points.

        Parameters
        ----------
        euler_grid_deg : (M, 3) Euler angles in degrees
        batch_size     : number of grid points per batch (memory control)

        Returns
        -------
        odf_values : (M,) ODF values in m.r.d.
        """
        from .kernel import de_la_vallee_poussin

        euler_rad = np.radians(euler_grid_deg)
        q_grid = _euler_to_quat(euler_rad[:, 0], euler_rad[:, 1], euler_rad[:, 2])
        q_grid = q_grid * np.sign(q_grid[:, 0:1] + 1e-30)

        M = len(q_grid)
        NS = len(self._q_sym_data)
        odf = np.empty(M, dtype=float)

        for start in range(0, M, batch_size):
            end = min(start + batch_size, M)
            q_b = q_grid[start:end]  # (B, 4)

            # cos(omega/2) = |dot(q_grid, q_sym_data)| for each (grid, data*sym) pair
            # (B, NS) = (B, 4) @ (4, NS)
            cos_half = np.abs(q_b @ self._q_sym_data.T)   # (B, NS)
            cos_half = np.clip(cos_half, 0.0, 1.0)

            # Reshape to (B, N, S) and SUM over sym ops, divide by S (MTEX convention)
            cos_half_3 = cos_half.reshape(len(q_b), self._N, self._S)
            omega_3 = 2.0 * np.arccos(cos_half_3)          # (B, N, S)
            k_3 = de_la_vallee_poussin(omega_3, self._kappa)  # (B, N, S)
            k_vals = k_3.mean(axis=2)                          # (B, N) — average over sym ops

            if self._weights is not None:
                odf[start:end] = (k_vals * self._weights[None, :]).sum(axis=1)
            else:
                odf[start:end] = k_vals.mean(axis=1)

        return odf

    # ------------------------------------------------------------------
    # ODF sections
    # ------------------------------------------------------------------

    def phi2_sections(self,
                      phi2_values_deg: list[float] | None = None,
                      resolution_deg: float = 5.0,
                      phi1_max_deg: float = 360.0,
                      PHI_max_deg: float = 90.0
                      ) -> dict[float, tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """
        Compute ODF for a set of φ₂ = const sections.

        Parameters
        ----------
        phi2_values_deg : list of φ₂ values in degrees; None → 0:5:90
        resolution_deg  : grid step for φ₁ and Φ
        phi1_max_deg    : maximum φ₁ value
        PHI_max_deg     : maximum Φ value

        Returns
        -------
        sections : dict mapping φ₂ → (phi1_2d, PHI_2d, odf_2d) all in degrees
                   phi1_2d and PHI_2d are meshgrid arrays; odf_2d is m.r.d.
        """
        if phi2_values_deg is None:
            phi2_values_deg = list(range(0, 91, 5))

        phi1_arr = np.arange(0, phi1_max_deg + 1e-9, resolution_deg)
        PHI_arr  = np.arange(0, PHI_max_deg  + 1e-9, resolution_deg)
        p1_2d, P_2d = np.meshgrid(phi1_arr, PHI_arr)

        sections = {}
        for phi2 in phi2_values_deg:
            p2_flat = np.full(p1_2d.size, phi2)
            grid = np.stack([p1_2d.ravel(), P_2d.ravel(), p2_flat], axis=1)
            odf_flat = self.evaluate(grid)
            odf_2d = odf_flat.reshape(p1_2d.shape)
            sections[float(phi2)] = (p1_2d, P_2d, odf_2d)

        return sections

    def sigma_sections(self,
                       sigma_values_deg: list[float] | None = None,
                       resolution_deg: float = 5.0
                       ) -> dict[float, tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """
        Sigma (φ₁ − φ₂ = const) sections.

        Returns
        -------
        sections : dict mapping sigma → (phi2_2d, PHI_2d, odf_2d)
        """
        if sigma_values_deg is None:
            sigma_values_deg = list(range(0, 91, 5))

        phi2_arr = np.arange(0, 91, resolution_deg)
        PHI_arr  = np.arange(0, 91, resolution_deg)
        p2_2d, P_2d = np.meshgrid(phi2_arr, PHI_arr)

        sections = {}
        for sigma in sigma_values_deg:
            p1_flat = (sigma + p2_2d.ravel()) % 360.0
            grid = np.stack([p1_flat, P_2d.ravel(), p2_2d.ravel()], axis=1)
            odf_flat = self.evaluate(grid)
            odf_2d = odf_flat.reshape(p2_2d.shape)
            sections[float(sigma)] = (p2_2d, P_2d, odf_2d)

        return sections

    def phi1_sections(self,
                      phi1_values_deg: list[float] | None = None,
                      resolution_deg: float = 5.0
                      ) -> dict[float, tuple[np.ndarray, np.ndarray, np.ndarray]]:
        """φ₁ = const sections."""
        if phi1_values_deg is None:
            phi1_values_deg = list(range(0, 91, 5))

        phi2_arr = np.arange(0, 91, resolution_deg)
        PHI_arr  = np.arange(0, 91, resolution_deg)
        p2_2d, P_2d = np.meshgrid(phi2_arr, PHI_arr)

        sections = {}
        for phi1 in phi1_values_deg:
            p1_flat = np.full(p2_2d.size, phi1)
            grid = np.stack([p1_flat, P_2d.ravel(), p2_2d.ravel()], axis=1)
            odf_flat = self.evaluate(grid)
            odf_2d = odf_flat.reshape(p2_2d.shape)
            sections[float(phi1)] = (p2_2d, P_2d, odf_2d)

        return sections

    # ------------------------------------------------------------------
    # Properties
    # ------------------------------------------------------------------

    def texture_index(self, resolution_deg: float = 5.0) -> float:
        """
        Texture index (J) = ∫ f(g)² dg  — 1 for random, >1 for textured.

        Estimated numerically on a grid.
        """
        phi1 = np.arange(0, 360, resolution_deg)
        PHI  = np.arange(0, 91,  resolution_deg)
        phi2 = np.arange(0, 360, resolution_deg)

        p1, P, p2 = np.meshgrid(phi1, PHI, phi2)
        grid = np.stack([p1.ravel(), P.ravel(), p2.ravel()], axis=1)

        odf = self.evaluate(grid, batch_size=500)
        P_rad = np.radians(P.ravel())
        weights = np.sin(P_rad)
        weights /= (weights.sum() + 1e-15)

        return float((odf ** 2 * weights).sum())

    def entropy(self, resolution_deg: float = 5.0) -> float:
        """
        Texture entropy = -∫ f(g) log(f(g)) dg
        Minimum for single-crystal, 0 for random.
        """
        phi1 = np.arange(0, 360, resolution_deg)
        PHI  = np.arange(0, 91,  resolution_deg)
        phi2 = np.arange(0, 360, resolution_deg)
        p1, P, p2 = np.meshgrid(phi1, PHI, phi2)
        grid = np.stack([p1.ravel(), P.ravel(), p2.ravel()], axis=1)
        odf  = self.evaluate(grid, batch_size=500)
        P_rad = np.radians(P.ravel())
        weights = np.sin(P_rad)
        weights /= (weights.sum() + 1e-15)
        f_pos = odf[odf > 0]
        w_pos = weights[odf > 0]
        return float(-(f_pos * np.log(f_pos) * w_pos).sum())

    def odf_max(self, resolution_deg: float = 5.0) -> float:
        """Maximum ODF value (m.r.d.)."""
        sections = self.phi2_sections(list(range(0, 91, resolution_deg)),
                                      resolution_deg)
        return max(v[2].max() for v in sections.values())

    # ------------------------------------------------------------------
    # Pole figures from ODF
    # ------------------------------------------------------------------

    def calc_pole_figure(self,
                         hkl: tuple = (0, 0, 1),
                         sample_symmetry: str = "none",
                         n_points: int = 2000
                         ) -> tuple[np.ndarray, np.ndarray | None]:
        """
        Compute pole figure vectors from the ODF (simulation mode).

        Samples random orientations from the ODF and projects crystal direction h.

        Returns
        -------
        vectors : (M, 3) unit vectors in sample frame (upper hemisphere)
        weights : None (uniform)
        """
        from .pole_figure import calc_pole_figure as _pf
        rand_euler = self.random_sampling(n_points)
        return _pf(rand_euler, hkl, self.crystal_symmetry, sample_symmetry)

    # ------------------------------------------------------------------
    # Random sampling from the ODF
    # ------------------------------------------------------------------

    def random_sampling(self, n: int, rng: np.random.Generator | None = None
                        ) -> np.ndarray:
        """
        Draw n orientations from the ODF using rejection sampling.

        Returns
        -------
        euler_deg : (n, 3) sampled Euler angles in degrees
        """
        if rng is None:
            rng = np.random.default_rng()

        # Pick from data orientations with probability ∝ ODF value at that point
        # Simple approach: pick data points (weighted) then perturb by kernel
        N = self._N
        if self._weights is not None:
            probs = self._weights / self._weights.sum()
        else:
            probs = np.ones(N) / N

        # Draw base orientations
        idx = rng.choice(N, size=n, p=probs)
        base_euler = self._euler_deg[idx]  # (n, 3)

        # Perturb by kernel halfwidth (Gaussian approximation)
        hw_rad = np.radians(self.halfwidth_deg) / 2.35  # FWHM → sigma
        noise = rng.normal(0, hw_rad, base_euler.shape)
        sampled = base_euler + np.degrees(noise)

        # Wrap phi1 and phi2 to [0, 360), clamp PHI to [0, 90]
        sampled[:, 0] = sampled[:, 0] % 360
        sampled[:, 1] = np.clip(sampled[:, 1], 0, 90)
        sampled[:, 2] = sampled[:, 2] % 360

        return sampled

    # ------------------------------------------------------------------
    # Volume fractions / component analysis
    # ------------------------------------------------------------------

    def volume_fraction(self, component_euler_deg: tuple,
                        max_angle_deg: float = 15.0) -> float:
        """Volume fraction within max_angle of a texture component."""
        from .components import volume_fraction as _vf
        return _vf(self._euler_deg, component_euler_deg, max_angle_deg,
                   self._weights)

    def component_analysis(self, max_angle_deg: float = 15.0) -> list[dict]:
        """Volume fractions for all standard components."""
        from .components import component_table
        return component_table(self._euler_deg, self.crystal_symmetry,
                                max_angle_deg, self._weights)

    # ------------------------------------------------------------------
    # Fibre ODF
    # ------------------------------------------------------------------

    def calc_fibre(self,
                   crystal_direction: tuple = (1, 0, 0),
                   sample_direction: tuple  = (0, 0, 1),
                   halfwidth_deg: float = 10.0) -> "FibreODF":
        """Return a FibreODF for the given crystal→sample direction pair."""
        return FibreODF(crystal_direction, sample_direction, halfwidth_deg,
                        self.crystal_symmetry)


# ---------------------------------------------------------------------------
# Fibre ODF
# ---------------------------------------------------------------------------

class FibreODF:
    """
    A model fibre ODF: all orientations g such that g·h is close to r.

    Parameters
    ----------
    h             : crystal direction (Miller indices, need not be normalised)
    r             : sample direction (need not be normalised)
    halfwidth_deg : angular FWHM of the fibre tube
    crystal_symmetry : crystal symmetry string
    """

    def __init__(self,
                 h: tuple = (1, 0, 0),
                 r: tuple = (0, 0, 1),
                 halfwidth_deg: float = 10.0,
                 crystal_symmetry: str = "orthorhombic"):
        self.h = np.array(h, dtype=float)
        self.h /= np.linalg.norm(self.h)
        self.r = np.array(r, dtype=float)
        self.r /= np.linalg.norm(self.r)
        self.halfwidth_deg = halfwidth_deg
        self.crystal_symmetry = crystal_symmetry

    def evaluate(self, euler_grid_deg: np.ndarray) -> np.ndarray:
        """Evaluate fibre ODF (m.r.d.) at Euler angle grid points."""
        from santex.grains.reconstruction import euler_to_rotation_matrix
        euler_rad = np.radians(euler_grid_deg)
        R = euler_to_rotation_matrix(euler_rad[:, 0], euler_rad[:, 1], euler_rad[:, 2])
        # g·h → sample frame
        g_h = np.einsum("nij,j->ni", R, self.h)   # (M, 3)
        # angle between g·h and r
        dot = np.clip(g_h @ self.r, -1, 1)
        omega = np.arccos(dot)
        # Gaussian kernel centred on the fibre
        sigma = np.radians(self.halfwidth_deg) / 2.35
        return np.exp(-0.5 * omega ** 2 / (sigma ** 2 + 1e-15))

    def sample(self, n: int = 1000, rng=None) -> np.ndarray:
        """Sample n orientations lying on the fibre."""
        if rng is None:
            rng = np.random.default_rng()
        sigma = np.radians(self.halfwidth_deg) / 2.35

        # Build one canonical orientation g0: h → r
        # Then rotate around r to fill the fibre
        euler_list = []
        for _ in range(n):
            # Random angle around the fibre axis r
            angle_around = rng.uniform(0, 2 * np.pi)
            # Perturbation from ideal fibre
            dev = rng.normal(0, sigma)
            # Build rotation from Rodrigues formula
            R0 = _align_vectors(self.h, self.r)
            Rz = _rot_around(self.r, angle_around)
            Rdev = _rot_around(
                _perp(self.r, rng),
                dev
            )
            R_total = Rdev @ Rz @ R0
            euler_list.append(_rot_to_euler_deg(R_total))

        return np.array(euler_list)


# ---------------------------------------------------------------------------
# Rotation geometry helpers for FibreODF.sample
# ---------------------------------------------------------------------------

def _align_vectors(v1, v2):
    """Return rotation matrix R such that R @ v1 = v2."""
    v1 = v1 / (np.linalg.norm(v1) + 1e-15)
    v2 = v2 / (np.linalg.norm(v2) + 1e-15)
    axis = np.cross(v1, v2)
    sin_a = np.linalg.norm(axis)
    cos_a = np.clip(np.dot(v1, v2), -1, 1)
    if sin_a < 1e-10:
        return np.eye(3) if cos_a > 0 else -np.eye(3)
    axis /= sin_a
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + sin_a * K + (1 - cos_a) * K @ K


def _rot_around(axis, angle):
    """Rotation matrix around a unit axis by angle (radians)."""
    axis = axis / (np.linalg.norm(axis) + 1e-15)
    c, s = np.cos(angle), np.sin(angle)
    K = np.array([[0, -axis[2], axis[1]],
                  [axis[2], 0, -axis[0]],
                  [-axis[1], axis[0], 0]])
    return np.eye(3) + s * K + (1 - c) * K @ K


def _perp(v, rng):
    """A random unit vector perpendicular to v."""
    r = rng.standard_normal(3)
    r -= r.dot(v) * v
    r /= (np.linalg.norm(r) + 1e-15)
    return r


def _rot_to_euler_deg(R):
    """Extract Bunge ZXZ Euler angles in degrees from a rotation matrix."""
    PHI = np.arccos(np.clip(R[2, 2], -1, 1))
    if np.abs(R[2, 2]) < 1 - 1e-10:
        phi1 = np.arctan2(R[0, 2], -R[1, 2]) % (2 * np.pi)
        phi2 = np.arctan2(R[2, 0],  R[2, 1]) % (2 * np.pi)
    else:
        phi1 = 0.0
        phi2 = np.arctan2(R[0, 1], R[0, 0]) % (2 * np.pi)
    return np.degrees([phi1, PHI, phi2])


# ---------------------------------------------------------------------------
# halfwidth → kappa (local alias to avoid circular import)
# ---------------------------------------------------------------------------

def _halfwidth_to_kappa(halfwidth_deg: float) -> float:
    from .kernel import halfwidth_to_kappa
    return halfwidth_to_kappa(halfwidth_deg)


# ---------------------------------------------------------------------------
# Convenience constructor from EBSD DataFrame
# ---------------------------------------------------------------------------

def odf_from_ebsd(data,
                  phase: int | None = None,
                  crystal_symmetry: str = "orthorhombic",
                  halfwidth_deg: float = 10.0,
                  max_orientations: int = 5000) -> ODF:
    """
    Build an ODF from an EBSD DataFrame.

    Parameters
    ----------
    data             : pd.DataFrame with Phase, Euler1, Euler2, Euler3 columns
    phase            : phase index to select; None = use all indexed pixels
    crystal_symmetry : crystal symmetry string
    halfwidth_deg    : kernel halfwidth
    max_orientations : cap on number of orientations used

    Returns
    -------
    odf : ODF object
    """
    if phase is not None:
        data = data[data["Phase"] == phase]
    else:
        data = data[data["Phase"] > 0]

    euler_deg = data[["Euler1", "Euler2", "Euler3"]].values
    return ODF(euler_deg, crystal_symmetry=crystal_symmetry,
               halfwidth_deg=halfwidth_deg, max_orientations=max_orientations)
