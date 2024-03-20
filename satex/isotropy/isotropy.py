#Adapted from Hacker et al 2003

import numpy as np
from scipy.optimize import root_scalar
import json
import os

class Isotropy:
    def __init__(self):
        try:
            json_file = 'materials.json'
            json_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', json_file)
            with open(json_path, 'r') as f:
                self.data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            raise FileNotFoundError(f"Error loading JSON file '{json_file}': {e}")

    def get_available_phases(self):
        try:
            for key in range(0,len(self.data)):
                print('######################Available Phases######################')
                print('Material id: ' + str(self.data[str(key)]['phase'] + '       ' + 'Material name: ' + str(self.data[str(key)]['name'])))
            return [self.data[key]['phase'] for key in self.data]
        except (KeyError, TypeError) as e:
            print(f"Error accessing data: {e}")
            return []

    def get_phase_constants(self, phase):
        try:
            phase_data = next(item for item in self.data.values() if item.get('phase') == phase or item.get('name') == phase)
            if phase_data:
                return {
                    'id': phase_data['phase'],
                    'name': phase_data['name'],
                    'rho0': phase_data['density_298K(kg/m3)'],
                    'ao': phase_data['expansivity(/K)_a0'],
                    'akt0': phase_data['isothermal_bulk_modulus(P)_KT'],
                    'dkdp': phase_data['dKt/dP(KT_prime)'],
                    'amu0': phase_data['shear_modulus(Pa)_G'],
                    'dmudp': phase_data['dG/dP(G_prime)'],
                    'gam': phase_data['dlnG/dlnr_gamma'],
                    'grun': phase_data['first_gruinesen_parameter_gth'],
                    'delt': phase_data['second_gruinesen_parameter_dT']
                }
            else:
                raise KeyError(f"Phase '{phase}' not found in the data.")
        except (StopIteration, KeyError) as e:
            print(f"Error accessing phase constants: {e}")
            return None

    def calculate_seismic_properties(self, phase, temperature, pressure, return_vp_vs_vbulk=False, return_aktout=False):
        """
        Calculate velocities at P,T for a mineral given parameters.
        
        Inputs:
        temperature (in K)
        pressure (in GPa)
        
        Outputs:
        
        density in kg/m3
        vp, vs, vbulk in km/s
        
        """
        
        if isinstance(temperature,(int,float)):
            if isinstance(pressure,(int,float)):
                method_calc = 'scalar'
            else:
                raise ValueError('While temperature is a scalar, pressure is an array.')
        else:
            if len(temperature) != len(pressure):
                raise ValueError('Length of the temperature and pressure arrays do not match!')
            else:
                method_calc = 'array'
        
        #Unit conversion from GPa
        
        pressure = pressure * 1e9 #conversion to GPa to Pa
        
        phase_constants = self.get_phase_constants(phase)

        if phase_constants:
            
            # Constants
            reference_temperature = 25.0
            absolute_temperature = 273.0
            tolerance = 1.e-8

            # input
            rho0 = phase_constants['rho0']
            ao = phase_constants['ao']
            akt0 = phase_constants['akt0']
            dkdp = phase_constants['dkdp']
            amu0 = phase_constants['amu0']
            dmudp = phase_constants['dmudp']
            gam = phase_constants['gam']
            grun = phase_constants['grun']
            delt = phase_constants['delt']

            # first, integrate with temperature at P=0
            alpht = ao * (1. - 10. / np.sqrt(absolute_temperature + temperature))
            phi = ao * (temperature - reference_temperature) - 20. * ao * (np.sqrt(temperature + absolute_temperature) - np.sqrt(reference_temperature + absolute_temperature))
            density_t = rho0 * np.exp(-phi)

            akt_t = akt0 * np.exp(-delt * phi)  # KtT in Excel macro
            amu_t = amu0 * np.exp(-gam * phi)   # GT in Excel macro

            # second, integrate with pressure at T fixed
            f_guess = pressure / (3. * akt_t)
            f_guess2 = 0.
            akk = akt_t
            akk_prime = dkdp

            if method_calc == 'array':

                ff = []
                
                for i in range(0,len(temperature)):

                    res = root_scalar(self.pressure_function, args=(pressure[i], akk[i], akk_prime), bracket=[f_guess[i], f_guess2], xtol=tolerance)
                    if not res.converged:
                        raise ValueError("Root finding did not converge at array index: " + str(i))
                    ff_ = res.root
                    ff.append(ff_)

                ff = np.array(ff)

            else:

                res = root_scalar(self.pressure_function, args=(pressure, akk, akk_prime), bracket=[f_guess, f_guess2], xtol=tolerance)
                if not res.converged:
                    raise ValueError("Root finding did not converge.")
                ff = res.root
            
            dvol = (1. + 2. * ff) ** 1.5
            density = density_t * dvol
            ffac = (1. + 2. * ff) ** 2.5

            akt = akt_t * ffac * (1. - (5. - 3. * dkdp) * ff)
            amu = amu_t * ffac * (1. - ff * (5. - 3. * dmudp * akt_t / amu_t))

            # convert to adiabatic bulk modulus
            alphtp = alpht * dvol ** (-delt)
            aks = akt * (1. + alphtp * grun * (temperature + absolute_temperature))
            
            if return_vp_vs_vbulk == True:
                vbulk, vp, vs = self.calculate_velocities(density, aks, amu)
                if return_aktout == True:
                    return density, aks, amu, vp, vs, vbulk, akt
                else:
                    return density, aks, amu, vp, vs, vbulk
            else:
                if return_aktout == True:
                    return density, aks, amu, akt
                else:
                    return density, aks, amu
        else:
            raise ValueError(f"No phase constants available for '{phase}'.")


    def calculate_velocities(self, density, aks, amu):
        try:
            vbulk = 0.001 * np.sqrt(aks / density)
            vp = 0.001 * np.sqrt((aks + 1.3333 * amu) / density)
            vs = 0.001 * np.sqrt(amu / density)
            return vbulk, vp, vs
        except (ZeroDivisionError, ValueError) as e:
            print(f"Error calculating velocities: {e}")

    def pressure_function(self, f, press, akk, akk_prime):
        try:
            zeta = 0.75 * (4. - akk_prime)
            f1 = (1. + 2. * f) ** 2.5
            f2 = (1. - 2. * zeta * f)
            pdif = 3. * akk * f * f1 * f2 - press
            return pdif
        except ZeroDivisionError as e:
            print(f"Error in pressure function: {e}")
                    
    def set_modal_composition(self, phase_list, fraction_list):
        
        """
        
        Inputs:
        
        phase_list - list of phase ids to be included in the rock. 
        Example: ['fo', 'en'] -- Forsterite and Enstatite
        
        fraction_list - list of fraction lists to be included in the rock
        Example: [0.8, 0.2] -- 0.8 fo and 0.2 en
        
        """
        
        phase_constant_list = []
        
        for item in phase_list:
        
            phase_constant = self.get_phase_constants(phase = item)
            phase_constant_list.append(phase_constant)
            
        return phase_constant_list, fraction_list
        
        