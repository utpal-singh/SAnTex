import numpy as np
from santex.anisotropy import Anisotropy
from santex.anisotropy.utils import christoffel_tensor, wave_property
from santex.tensor import Tensor


class AnisotropyBackend:
    def __init__(self):
        self.anisotropy: Anisotropy | None = None
        self.cijkl: np.ndarray | None = None
        self.rho: float | None = None

    def set_from_voigt(self, voigt_pa: np.ndarray, density: float) -> None:
        """Initialise from a 6×6 Voigt matrix (Pa) and density (kg/m³)."""
        self.anisotropy = Anisotropy(voigt_pa, density)
        self.cijkl = self.anisotropy.cijkl
        self.rho = density

    # ------------------------------------------------------------------
    # Scalar anisotropy metrics
    # ------------------------------------------------------------------

    def compute_anisotropy_values(self) -> dict | None:
        if self.anisotropy is None:
            return None
        return self.anisotropy.anisotropy_values()

    def voigt_velocity(self):
        return self.anisotropy.voigt_velocity() if self.anisotropy else None

    def reuss_velocity(self):
        return self.anisotropy.reuss_velocity() if self.anisotropy else None

    def hill_velocity(self):
        return self.anisotropy.hill_velocity() if self.anisotropy else None

    # ------------------------------------------------------------------
    # Vectorised velocity surface (fast – uses einsum, no Python loops)
    # ------------------------------------------------------------------

    def compute_velocity_surface(self, n_theta: int = 60, n_phi: int = 120) -> dict | None:
        """Return velocity arrays over a (n_theta × n_phi) spherical grid."""
        if self.cijkl is None or self.rho is None:
            return None

        theta = np.linspace(0, np.pi, n_theta)
        phi = np.linspace(0, 2 * np.pi, n_phi)
        th, ph = np.meshgrid(theta, phi, indexing="ij")

        nx = np.sin(th) * np.cos(ph)
        ny = np.sin(th) * np.sin(ph)
        nz = np.cos(th)

        # Flatten directions → (N, 3)
        n_flat = np.stack([nx.ravel(), ny.ravel(), nz.ravel()], axis=1)

        # Vectorised Christoffel tensor: T_ik = C_ijkl n_j n_l → shape (N, 3, 3)
        tik = np.einsum("ijkl,nj,nl->nik", self.cijkl, n_flat, n_flat)

        # eigh guarantees real eigenvalues for symmetric matrices; returns ascending
        evals = np.linalg.eigh(tik)[0][:, ::-1]  # (N, 3) descending

        vp  = np.sqrt(np.maximum(evals[:, 0], 0) / self.rho).reshape(n_theta, n_phi)
        vs1 = np.sqrt(np.maximum(evals[:, 1], 0) / self.rho).reshape(n_theta, n_phi)
        vs2 = np.sqrt(np.maximum(evals[:, 2], 0) / self.rho).reshape(n_theta, n_phi)

        denom_avs  = vs1 + vs2
        denom_vvs1 = vs1
        denom_vvs2 = vs2
        avs   = np.where(denom_avs  > 0, 200 * (vs1 - vs2) / denom_avs,  0.0)
        vpvs1 = np.where(denom_vvs1 > 0, vp / denom_vvs1, 0.0)
        vpvs2 = np.where(denom_vvs2 > 0, vp / denom_vvs2, 0.0)

        return {
            "nx": nx, "ny": ny, "nz": nz,
            "vp": vp, "vs1": vs1, "vs2": vs2,
            "avs": avs, "vpvs1": vpvs1, "vpvs2": vpvs2,
        }

    # ------------------------------------------------------------------
    # Upper-hemisphere stereonet data for 2-D plots
    # ------------------------------------------------------------------

    def compute_stereonet_data(self, step_deg: float = 3.0) -> dict | None:
        """Return scatter data in stereographic projection (upper hemisphere)."""
        if self.cijkl is None or self.rho is None:
            return None

        step = np.deg2rad(step_deg)
        theta_vals = np.arange(0, np.pi / 2 + step, step)
        phi_vals   = np.arange(0, 2 * np.pi + step, step)
        th, ph = np.meshgrid(theta_vals, phi_vals, indexing="ij")

        nx = np.sin(th) * np.cos(ph)
        ny = np.sin(th) * np.sin(ph)
        nz = np.cos(th)
        n_flat = np.stack([nx.ravel(), ny.ravel(), nz.ravel()], axis=1)

        tik   = np.einsum("ijkl,nj,nl->nik", self.cijkl, n_flat, n_flat)
        evals = np.linalg.eigh(tik)[0][:, ::-1]

        vp  = np.sqrt(np.maximum(evals[:, 0], 0) / self.rho)
        vs1 = np.sqrt(np.maximum(evals[:, 1], 0) / self.rho)
        vs2 = np.sqrt(np.maximum(evals[:, 2], 0) / self.rho)

        # Stereographic projection: x = n_x/(1+n_z), y = n_y/(1+n_z)
        nz_flat = nz.ravel()
        denom   = 1 + nz_flat
        x = nx.ravel() / denom
        y = ny.ravel() / denom

        avs   = np.where((vs1 + vs2) > 0, 200 * (vs1 - vs2) / (vs1 + vs2), 0.0)
        vpvs1 = np.where(vs1 > 0, vp / vs1, 0.0)
        vpvs2 = np.where(vs2 > 0, vp / vs2, 0.0)

        return {
            "x": x, "y": y,
            "vp": vp, "vs1": vs1, "vs2": vs2,
            "avs": avs, "vpvs1": vpvs1, "vpvs2": vpvs2,
        }
