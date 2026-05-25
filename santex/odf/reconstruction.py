"""
ODF Reconstruction from pole figures (iterative WIMV-like algorithm).

The WIMV (Williams–Imhof–Matthies–Vinel) algorithm reconstructs f(g) from
measured pole figures P_h(y) by iterating:

  f^(n+1)(g) ← f^(n)(g) · ∏_h [ P_h(g·h) / P_h^(n)(g·h) ]^(1/|h|)

where P_h^(n) is the recalculated pole figure from f^(n).

Reference: Matthies, Vinel & Helming (1987), "Standard Distributions in
Texture Analysis", Vol. 1.
"""

from __future__ import annotations
import numpy as np


class PoleFigureReconstruction:
    """
    Iterative ODF reconstruction from experimental pole figures.

    Parameters
    ----------
    crystal_symmetry : crystal symmetry string
    resolution_deg   : grid resolution for the ODF (5° recommended)
    n_iter           : number of WIMV iterations
    """

    def __init__(self,
                 crystal_symmetry: str = "cubic",
                 resolution_deg: float = 5.0,
                 n_iter: int = 20):
        self.crystal_symmetry = crystal_symmetry
        self.resolution_deg   = resolution_deg
        self.n_iter           = n_iter
        self._pole_figures: list[tuple] = []   # (hkl, pf_vectors, pf_intensities)

    # ------------------------------------------------------------------

    def add_pole_figure(self, hkl: tuple, vectors: np.ndarray,
                        intensities: np.ndarray) -> None:
        """
        Register an experimental pole figure.

        Parameters
        ----------
        hkl         : Miller indices of the reflection, e.g. (1, 0, 0)
        vectors     : (N, 3) measurement directions in sample frame (unit vectors)
        intensities : (N,) measured intensity values (normalised to m.r.d.)
        """
        h = np.array(hkl, dtype=float)
        h /= np.linalg.norm(h) + 1e-15
        self._pole_figures.append((h, np.asarray(vectors, float),
                                   np.asarray(intensities, float)))

    def reconstruct(self) -> np.ndarray:
        """
        Run WIMV iterations and return ODF values on the internal grid.

        Returns
        -------
        odf : (M,) ODF values in m.r.d. at the grid points (phi1, PHI, phi2)
        grid : (M, 3) corresponding Euler angle grid in degrees
        """
        from santex.grains.reconstruction import euler_to_rotation_matrix
        from santex.grains.symmetry import SYMMETRY

        res = self.resolution_deg
        phi1 = np.arange(0, 360, res)
        PHI  = np.arange(0, 91,  res)
        phi2 = np.arange(0, 360, res)
        p1g, Pg, p2g = np.meshgrid(phi1, PHI, phi2)
        grid = np.stack([p1g.ravel(), Pg.ravel(), p2g.ravel()], axis=1)  # (M, 3)
        M = len(grid)

        euler_rad = np.radians(grid)
        R = euler_to_rotation_matrix(
            euler_rad[:, 0], euler_rad[:, 1], euler_rad[:, 2]
        )  # (M, 3, 3)

        csym = SYMMETRY.get(self.crystal_symmetry.lower(), SYMMETRY["triclinic"])

        # Precompute rotation of h through each R and each crystal sym op
        # For each pf: pf_dir[m, sc] = R[m] @ csym[sc] @ h  → (M, S, 3)

        # Initialise ODF to uniform (1 m.r.d.)
        f = np.ones(M, dtype=float)

        for iteration in range(self.n_iter):
            f_new = np.ones(M, dtype=float)

            for h_vec, y_vecs, I_exp in self._pole_figures:
                # Recalculate pole figure from current ODF
                # P_h(y) = (1/|cs|) Σ_{sc} ∫ f(g) δ(g·h_sc - y) dg
                # Discretised: for each measurement direction y_i,
                # find grid orientations g where g·h ≈ y_i

                # Apply crystal symmetry to h
                h_sym = (csym @ h_vec)    # (S, 3)
                S = len(h_sym)

                # g·h for all grid orientations and all sym-equivalent h:
                # g_h[m, sc] = R[m] @ h_sym[sc]
                g_h = np.einsum("mij,sj->msi", R, h_sym)  # (M, S, 3)
                g_h_flat = g_h.reshape(-1, 3)              # (M*S, 3)

                # For each measurement direction y_i, accumulate intensity
                # from grid points where the angle g·h_sc ≈ y_i (within res/2)
                I_calc = np.zeros(len(y_vecs))
                for yi, y in enumerate(y_vecs):
                    # Dot product with measurement directions
                    dots = np.clip(g_h_flat @ y, -1, 1)
                    omega = np.arccos(dots)  # angle to measurement dir
                    within = omega < np.radians(res)
                    if within.any():
                        # weights from f values
                        m_idx = np.where(within)[0] // S
                        I_calc[yi] = f[m_idx].mean()
                    else:
                        I_calc[yi] = 1.0

                # WIMV update ratio: P_exp / P_calc for each orientation
                for m in range(M):
                    # Find which measurement directions are close to R[m]@h
                    for sc in range(S):
                        d = R[m] @ h_sym[sc]
                        dots = np.clip(y_vecs @ d, -1, 1)
                        omega = np.arccos(dots)
                        i_near = np.argmin(omega)
                        if I_calc[i_near] > 0:
                            f_new[m] *= (I_exp[i_near] / I_calc[i_near]) ** (1.0 / len(self._pole_figures))

            f = f * f_new
            # Normalise
            P_rad = np.radians(grid[:, 1])
            weights = np.sin(P_rad)
            f /= (f * weights).mean() + 1e-15

        return f, grid

    def ghost_correction(self, odf_values: np.ndarray,
                         method: str = "positivity") -> np.ndarray:
        """
        Apply ghost correction to an ODF estimated from pole figures.

        The ghost effect causes negative ODF values (odd-rank harmonics are
        zero from pole figures) which are then set to zero, creating the
        'ghost'. Positivity correction: set negative values to zero and
        renormalise.

        Parameters
        ----------
        odf_values : (M,) ODF values from reconstruction
        method     : 'positivity' (PhilSCF method) — only option currently

        Returns
        -------
        corrected : (M,) ghost-corrected ODF values
        """
        corrected = np.maximum(odf_values, 0.0)
        mean_val = corrected.mean()
        if mean_val > 0:
            corrected /= mean_val
        return corrected
