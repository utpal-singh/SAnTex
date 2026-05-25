import numpy as np
from santex.material import Material
from santex.isotropy import Isotropy


def _to_float(val, default=0.0) -> float:
    try:
        return float(val)
    except (TypeError, ValueError):
        return default


class MaterialBackend:
    def __init__(self):
        self.material = Material()
        try:
            self.isotropy = Isotropy()
        except Exception:
            self.isotropy = None

    # ------------------------------------------------------------------
    # Material list helpers
    # ------------------------------------------------------------------

    def get_all_phases(self) -> list[dict]:
        """Return list of dicts with keys: phase, crystal_system, primary_phase, density."""
        out = []
        for m in self.material.materials_data:
            out.append({
                "phase":          m.get("Phase", ""),
                "crystal_system": m.get("Crystal System", ""),
                "primary_phase":  m.get("Primary Phase", ""),
                "density":        _to_float(m.get("Density(g/cm3)", 0.0)),
            })
        return out

    def get_phase_names(self) -> list[str]:
        return [m.get("Phase", "") for m in self.material.materials_data]

    # ------------------------------------------------------------------
    # Stiffness / density at P, T
    # ------------------------------------------------------------------

    def get_voigt_matrix_gpa(self, phase: str, pressure_gpa: float = 0.0,
                              temp_k: float = 300.0) -> np.ndarray | None:
        """Return 6×6 Voigt stiffness matrix in GPa."""
        try:
            return self.material.voigt_high_PT(phase, PRESSURE=pressure_gpa, TEMP=temp_k)
        except Exception as e:
            print(f"MaterialBackend.get_voigt_matrix_gpa: {e}")
            return None

    def get_density(self, phase: str) -> float | None:
        """Return density in kg/m³ (from database, ambient conditions)."""
        try:
            return self.material.load_density(phase)
        except Exception as e:
            print(f"MaterialBackend.get_density: {e}")
            return None

    def get_raw_voigt_gpa(self, phase: str) -> np.ndarray | None:
        """Return the ambient-conditions 6×6 Voigt matrix in GPa."""
        try:
            return self.material.get_voigt_matrix(phase)
        except Exception as e:
            print(f"MaterialBackend.get_raw_voigt_gpa: {e}")
            return None

    # ------------------------------------------------------------------
    # Modal / multi-mineral rock
    # ------------------------------------------------------------------

    def modal_rock(self, minerals: list[str], fractions: list[float],
                   pressure_gpa: float, temp_k: float) -> tuple[np.ndarray | None, float | None]:
        """Voigt average for a multi-mineral assemblage. Returns (cij_GPa, density_kg_m3)."""
        try:
            cij, rho = self.material.modal_rock(minerals, fractions, pressure_gpa, temp_k)
            return cij, rho
        except Exception as e:
            print(f"MaterialBackend.modal_rock: {e}")
            return None, None

    # ------------------------------------------------------------------
    # Isotropic properties via Birch-Murnaghan EOS
    # ------------------------------------------------------------------

    def get_isotropic_phases(self) -> list[tuple[str, str]]:
        """Return list of (phase_id, phase_name) available in the isotropy database."""
        if self.isotropy is None:
            return []
        try:
            return [(v.get("phase", k), v.get("name", k))
                    for k, v in self.isotropy.data.items()]
        except Exception:
            return []

    def calc_isotropic_velocities(self, phase_id: str, temp_c: float,
                                  pressure_gpa: float) -> dict | None:
        """
        Returns dict with density, vp, vs, vbulk (all in SI / km·s⁻¹).
        temperature in °C, pressure in GPa.
        """
        if self.isotropy is None:
            return None
        try:
            density, aks, amu, vp, vs, vbulk = self.isotropy.calculate_seismic_properties(
                phase_id, temperature=temp_c, pressure=pressure_gpa, return_vp_vs_vbulk=True
            )
            return {"density": density, "vp": vp, "vs": vs, "vbulk": vbulk,
                    "K_adiabatic": aks, "G": amu}
        except Exception as e:
            print(f"MaterialBackend.calc_isotropic_velocities: {e}")
            return None

    def hashin_shtrikman(self, phase_ids: list[str], fractions: list[float],
                          temp_c: float, pressure_gpa: float) -> dict | None:
        """Return HS upper, lower and VRH bounds."""
        if self.isotropy is None:
            return None
        try:
            phase_constants, frac_arr = self.isotropy.set_modal_composition(phase_ids, fractions)
            medium, upper, lower, rho_list, K_vrh, G_vrh = self.isotropy.hashin_shtrikman_bounds(
                phase_constants, frac_arr, temp_c, pressure_gpa,
                density_mix_calc=True, modulii_return=True,
            )
            return {"medium": medium, "upper": upper, "lower": lower,
                    "density": rho_list, "K_vrh": K_vrh, "G_vrh": G_vrh}
        except Exception as e:
            print(f"MaterialBackend.hashin_shtrikman: {e}")
            return None
