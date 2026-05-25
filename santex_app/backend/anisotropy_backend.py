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
    # Upper-hemisphere stereonet — regular Cartesian grid (for Heatmap)
    # ------------------------------------------------------------------

    def compute_stereonet_grid(
        self, grid_size: int = 300, step_deg: float = 1.0
    ) -> dict | None:
        """
        Return velocity data on a ``grid_size × grid_size`` Cartesian grid
        using the **Lambert equal-area** (Schmidt) projection — the same
        default used by MTEX.  Values outside the unit disk are NaN.

        Projection:  r = √2 · sin(θ/2),  inverse: θ = 2·arcsin(r/√2)
        This preserves area (unlike the equal-angle/Wulff projection) so
        the visual weight of each solid-angle element is constant.

        Uses bilinear interpolation on the regular (theta, phi) source grid
        — no scattered-data interpolation, no edge artefacts.
        """
        if self.cijkl is None or self.rho is None:
            return None

        from scipy.ndimage import map_coordinates

        step = np.deg2rad(step_deg)
        theta_arr = np.arange(0.0, np.pi / 2 + step, step)
        phi_arr   = np.arange(0.0, 2 * np.pi + step, step)
        n_th, n_ph = len(theta_arr), len(phi_arr)

        th, ph = np.meshgrid(theta_arr, phi_arr, indexing="ij")
        nx_g = np.sin(th) * np.cos(ph)
        ny_g = np.sin(th) * np.sin(ph)
        nz_g = np.cos(th)
        n_flat = np.stack([nx_g.ravel(), ny_g.ravel(), nz_g.ravel()], axis=1)

        tik   = np.einsum("ijkl,nj,nl->nik", self.cijkl, n_flat, n_flat)
        evals = np.linalg.eigh(tik)[0][:, ::-1]   # descending → Vp, Vs1, Vs2

        vp  = np.sqrt(np.maximum(evals[:, 0], 0) / self.rho).reshape(n_th, n_ph)
        vs1 = np.sqrt(np.maximum(evals[:, 1], 0) / self.rho).reshape(n_th, n_ph)
        vs2 = np.sqrt(np.maximum(evals[:, 2], 0) / self.rho).reshape(n_th, n_ph)
        avs   = np.where((vs1 + vs2) > 0, 200 * (vs1 - vs2) / (vs1 + vs2), 0.0)
        vpvs1 = np.where(vs1 > 0, vp / vs1, 0.0)
        vpvs2 = np.where(vs2 > 0, vp / vs2, 0.0)

        # ── Lambert equal-area output grid ──────────────────────────────
        # Disk radius r ∈ [0,1] ↔ θ ∈ [0, 90°]
        # Forward:  r = √2 · sin(θ/2)
        # Inverse:  θ = 2 · arcsin(r / √2)   (valid for r ≤ 1 since r/√2 ≤ 1/√2 < 1)
        xi = np.linspace(-1.0, 1.0, grid_size)
        yi = np.linspace(-1.0, 1.0, grid_size)
        XI, YI = np.meshgrid(xi, yi)                    # (grid_size, grid_size)

        r2      = XI ** 2 + YI ** 2
        outside = r2 > 1.0
        r_safe  = np.sqrt(np.where(outside, 0.0, r2))

        # Inverse Lambert equal-area: r = √2·sin(θ/2)  →  θ = 2·arcsin(r/√2)
        theta_out = 2.0 * np.arcsin(np.clip(r_safe / np.sqrt(2.0), 0.0, 1.0))
        phi_out   = np.arctan2(YI, XI) % (2 * np.pi)   # [0, 2π)

        # Convert to fractional source-grid indices
        th_idx = theta_out / (np.pi / 2) * (n_th - 1)
        ph_idx = phi_out   / (2 * np.pi) * (n_ph - 1)

        coords = [th_idx.ravel(), ph_idx.ravel()]

        result = {"xi": xi, "yi": yi}
        for key, src in [("vp", vp), ("vs1", vs1), ("vs2", vs2),
                         ("avs", avs), ("vpvs1", vpvs1), ("vpvs2", vpvs2)]:
            grid = map_coordinates(src, coords, order=2,
                                   mode="nearest").reshape(grid_size, grid_size)
            grid[outside] = np.nan
            result[key] = grid

        return result

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

        # Lambert equal-area projection (same as MTEX default):
        #   x = √(2/(1+nz)) · nx,  y = √(2/(1+nz)) · ny
        # This is equivalent to  r = √2·sin(θ/2),  x = r·cos(φ),  y = r·sin(φ)
        nz_flat = nz.ravel()
        scale   = np.sqrt(np.maximum(2.0 / (1.0 + nz_flat), 0.0))
        x = scale * nx.ravel()
        y = scale * ny.ravel()

        avs   = np.where((vs1 + vs2) > 0, 200 * (vs1 - vs2) / (vs1 + vs2), 0.0)
        vpvs1 = np.where(vs1 > 0, vp / vs1, 0.0)
        vpvs2 = np.where(vs2 > 0, vp / vs2, 0.0)

        return {
            "x": x, "y": y,
            "vp": vp, "vs1": vs1, "vs2": vs2,
            "avs": avs, "vpvs1": vpvs1, "vpvs2": vpvs2,
        }
