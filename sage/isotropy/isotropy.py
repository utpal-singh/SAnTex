import numpy as np
from scipy.optimize import root_scalar

class Isotropy:
    def __init__(self, mineral_properties):
        self.mineral_properties = mineral_properties

    def calculate_seismic_properties(self, temperature, pressure, return_vp_vs_vbulk=False, return_aktout=False):
        """
        Calculate velocities at P,T for a mineral given parameters.
        """
        # Special handling if this is quartz
        if self.mineral_properties['name'] in ['aqz', 'bqz', 'qz']:
            density, aks, amu, akt = ah16_mineral_quartz(temperature, pressure, self.mineral_properties)
            if return_vp_vs_vbulk:
                vbulk, vp, vs = self.calculate_velocities(density, aks, amu)
                if return_aktout:
                    return density, aks, amu, vp, vs, vbulk, akt
                return density, aks, amu, vp, vs, vbulk
            return density, aks, amu

        # Constants
        reference_temperature = 25.0
        absolute_temperature = 273.0
        tolerance = 1.e-8

        # input
        rho0 = self.mineral_properties['rho0']
        ao = self.mineral_properties['ao']
        akt0 = self.mineral_properties['akt0']
        dkdp = self.mineral_properties['dkdp']
        amu0 = self.mineral_properties['amu0']
        dmudp = self.mineral_properties['dmudp']
        gam = self.mineral_properties['gam']
        grun = self.mineral_properties['grun']
        delt = self.mineral_properties['delt']

        # first, integrate with temperature at P=0
        alpht = ao * (1. - 10. / np.sqrt(absolute_temperature + temperature))
        phi = ao * (temperature - reference_temperature) - 20. * ao * (np.sqrt(temperature + absolute_temperature) - np.sqrt(reference_temperature + absolute_temperature))
        density_t = rho0 * np.exp(-phi)

        akt_t = akt0 * np.exp(-delt * phi)  # KtT in Excel macro
        amu_t = amu0 * np.exp(-gam * phi)   # GT in Excel macro

        # second, integrate with pressure at T fixed
        f_guess = pressure / (3. * akt_t)
        f_guess2 = 0.
        press = pressure
        akk = akt_t
        akk_prime = dkdp

        res = root_scalar(self.pressure_function, args=(press, akk, akk_prime), bracket=[f_guess, f_guess2], xtol=tolerance)
        if not res.converged:
            raise ValueError("Root finding did not converge")
        
        ff = res.root
        dvol = (1. + 2. * ff) ** 1.5
        density = density_t * dvol
        ffac = (1. + 2. * ff) ** 2.5

        akt = akt_t * ffac * (1. - (5. - 3. * dkdp) * ff)
        amu = amu_t * ffac * (1. - ff * (5. - 3. * dmudp * akt_t / amu_t))

        # convert to adiabatic bulk modulus
        alphtp = alpht * dvol ** (-delt)
        aks = akt * (1. + alphtp * grun * (temperature + absolute_temperature))

        if return_vp_vs_vbulk:
            vbulk, vp, vs = self.calculate_velocities(density, aks, amu)
            if return_aktout:
                return density, aks, amu, vp, vs, vbulk, akt
            return density, aks, amu, vp, vs, vbulk

        return density, aks, amu

    def calculate_velocities(self, density, aks, amu):
        vbulk = 0.001 * np.sqrt(aks / density)
        vp = 0.001 * np.sqrt((aks + 1.3333 * amu) / density)
        vs = 0.001 * np.sqrt(amu / density)
        return vbulk, vp, vs

    def pressure_function(self, f, press, akk, akk_prime):
        """
        Pressure vs strain function used for minimization (root finding).
        """
        zeta = 0.75 * (4. - akk_prime)
        f1 = (1. + 2. * f) ** 2.5
        f2 = (1. - 2. * zeta * f)
        pdif = 3. * akk * f * f1 * f2 - press
        return pdif

if __name__ == "__main__":
    mineral_properties = {
        'name': 'coe',
        'rho0': 2911,
        'ao': 1.8e-5,
        'akt0': 9.74e10,
        'dkdp': 4.3,
        'amu0': 6.16e10,
        'dmudp': 1.0541,
        'gam': 4.66,
        'grun': 0.36,
        'delt': 4.66,
    }

    isotropy = Isotropy(mineral_properties)
    density, aks, amu, vp, vs, vbulk, akt = isotropy.calculate_seismic_properties(2000, 2e9, return_vp_vs_vbulk=True, return_aktout=True)
    print(vp, vs)
