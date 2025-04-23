#Adapted from Hacker et al 2003

import numpy as np
from scipy.optimize import root_scalar
import json
import os

class Isotropy:
    def __init__(self):
        """
        Initializes the Isotropy class by loading material data from a JSON file.

        Raises:
        FileNotFoundError: If the JSON file 'materials.json' is not found.
        json.JSONDecodeError: If there is an error decoding the JSON file.
        """
        try:
            json_file = 'materials.json'
            json_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data', json_file)
            with open(json_path, 'r') as f:
                self.data = json.load(f)
        except (FileNotFoundError, json.JSONDecodeError) as e:
            raise FileNotFoundError(f"Error loading JSON file '{json_file}': {e}")

    def get_available_phases(self):
        """
        Prints and returns a list of available phases from the loaded material data.

        Returns:
        List: A list of available phases.
        """
        try:
            for key in range(0,len(self.data)):
                print('######################Available Phases######################')
                print('Material id: ' + str(self.data[str(key)]['phase'] + '       ' + 'Material name: ' + str(self.data[str(key)]['name'])))
            return [self.data[key]['phase'] for key in self.data]
        except (KeyError, TypeError) as e:
            print(f"Error accessing data: {e}")
            return []

    def get_phase_constants(self, phase):
        """
        Retrieves phase constants for a given phase name or ID from the loaded material data.

        Parameters:
        phase (str): The name or ID of the phase.

        Returns:
        dict: Phase constants including ID, name, density, expansivity, bulk modulus, shear modulus, and more.
        """
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

    def calculate_seismic_properties(self, phase, temperature, pressure, ref_density = None, return_vp_vs_vbulk=False, return_aktout=False):
        """
        Calculate velocities at P,T for a mineral given parameters.
        
        Inputs:
        temperature (in C)
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
            if ref_density == None:
                rho0 = phase_constants['rho0']
            else:
                rho0 = ref_density
            
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
        """
        Calculates velocities (vp, vs, vbulk) based on density, adiabatic bulk modulus (aks), and shear modulus (amu).

        Parameters:
        density (float or array): Density in kg/m^3.
        aks (float or array): Adiabatic bulk modulus in Pa.
        amu (float or array): Shear modulus in Pa.

        Returns:
        Tuple: Velocities vp, vs, vbulk in km/s.
        """
        try:
            vbulk = 0.001 * np.sqrt(aks / density)
            vp = 0.001 * np.sqrt((aks + 1.3333 * amu) / density)
            vs = 0.001 * np.sqrt(amu / density)
            return vbulk, vp, vs
        except (ZeroDivisionError, ValueError) as e:
            print(f"Error calculating velocities: {e}")

    def pressure_function(self, f, press, akk, akk_prime):
        """
        Helper function for pressure calculation in calculate_seismic_properties method.

        Parameters:
        f (float): Pressure factor.
        press (float): Pressure in Pa.
        akk (float): Adiabatic bulk modulus in Pa.
        akk_prime (float): Derivative of bulk modulus with respect to pressure.

        Returns:
        float: Pressure difference.
        """
        try:
            zeta = 0.75 * (4. - akk_prime)
            f1 = (1. + 2. * f) ** 2.5
            f2 = (1. - 2. * zeta * f)
            pdif = 3. * akk * f * f1 * f2 - press
            return pdif
        except ZeroDivisionError as e:
            print(f"Error in pressure function: {e}")
            
    def hashin_shtrikman_bounds(self, phase_constant_list, fraction_list, temperature, pressure, fractions_matches_T = False, density_mix_calc = False, modulii_return = False):
    
        """
        Calculates Hashin-Shtrikman Bounds for the given parameters:
        
        Parameters:
        - phase_constant_list: list of phase constant instances gotten from set_modal_composition function.
        - fraction_list: list of fractions of each phase.
        - temperature: temperature array or scalar in Kelvin
        - pressure: pressure array or scalar in GPa
        
        Returns:
        - 3 lists of medium value, upper bound and lower bound, indexed respectively: 
            [Vbulk_medium,Vp_medium,Vs_medium], [Vbulk_upper,Vp_upper,Vs_upper], [Vbulk_lower, Vp_lower, Vs_lower]
        
        """
        aks_list = []
        amu_list = []
        arho_list = []
        
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
                
                if len(fraction_list) == len(temperature):
                    fractions_matches_T = True
                else:
                    fractions_matches_T = False
                
        if method_calc == 'array':
            density_mix = np.zeros(len(temperature))
        else:
            density_mix = 0.0
        
        for phase_idx in range(0,len(phase_constant_list)):
        
            #calculating each phase moduli and density at T and P
            density, aks, amu = self.calculate_seismic_properties(phase_constant_list[phase_idx]['id'], temperature=temperature, pressure=pressure,
            return_vp_vs_vbulk=False, return_aktout=False)
            
            #calculating density of the mixture
            if fractions_matches_T == False:
                density_mix = density_mix + (density * fraction_list[phase_idx])
            else:
                density_mix = density_mix + (density * fraction_list[:,phase_idx])
            
            #appending moduli into lists
            aks_list.append(aks)
            amu_list.append(amu)
            
            if density_mix_calc == True:
                arho_list.append(density_mix)
            
        
        aks_list = np.array(aks_list)
        amu_list = np.array(amu_list)
        
        #Minimum maximum bounds calculation for moduli
        aks_min_array = np.amin(aks_list, axis = 0)
        aks_max_array = np.amax(aks_list, axis = 0)
        
        amu_min_array = np.amin(amu_list, axis = 0)
        amu_max_array = np.amax(amu_list, axis = 0)
        
        zet_min = (amu_min_array / 6.0) * ((9.0*aks_min_array) + (8.0*amu_min_array)) / (aks_min_array + (2*amu_min_array))
        zet_max = (amu_max_array / 6.0) * ((9.0*aks_max_array) + (8.0*amu_max_array)) / (aks_max_array + (2*amu_max_array))
        
        if method_calc == 'array':
            skm = np.zeros(len(temperature))
            skp = np.zeros(len(temperature))
            smm = np.zeros(len(temperature))
            smp = np.zeros(len(temperature))
        else:
            skm = 0.0
            skp = 0.0
            smm = 0.0
            smp = 0.0
        
        for phase_idx in range(0,len(phase_constant_list)):
            if fractions_matches_T == False:
                frac_list_holder = fraction_list[phase_idx]
            else:
                frac_list_holder = fraction_list[:,phase_idx]
                
            skm = skm + (frac_list_holder / (aks_list[phase_idx] + (4*amu_min_array/3)))
            skp = skp + (frac_list_holder / (aks_list[phase_idx] + (4*amu_max_array/3)))
            smm = smm + (frac_list_holder / (amu_list[phase_idx] + zet_min))
            smp = smp + (frac_list_holder / (amu_list[phase_idx] + zet_max))
            
        khsm = (1.0 / skm) - (4.0*(amu_min_array/3.0))
        khsp = (1.0 / skp) - (4.0*(amu_max_array/3.0))
        mhsm = (1.0/smm) - zet_min
        mhsp = (1.0/smp) - zet_max
        
        rkkh = 0.5*(khsm+khsp)
        rkgh = 0.5*(mhsm+mhsp)
        
        Vbulk_medium, Vp_medium, Vs_medium = self.calculate_velocities(density = density_mix, aks = rkkh ,amu = rkgh)
        Vbulk_upper, Vp_upper, Vs_upper = self.calculate_velocities(density = density_mix, aks = khsp ,amu = mhsp)
        Vbulk_lower, Vp_lower, Vs_lower = self.calculate_velocities(density = density_mix, aks = khsm ,amu = mhsm)
        
        if density_mix_calc == True:
            if modulii_return == True:
                return [Vbulk_medium,Vp_medium,Vs_medium], [Vbulk_upper,Vp_upper,Vs_upper], [Vbulk_lower, Vp_lower, Vs_lower], arho_list, rkkh, rkgh
            else:
                return [Vbulk_medium,Vp_medium,Vs_medium], [Vbulk_upper,Vp_upper,Vs_upper], [Vbulk_lower, Vp_lower, Vs_lower], arho_list
        else:
            if modulii_return == True:
                return [Vbulk_medium,Vp_medium,Vs_medium], [Vbulk_upper,Vp_upper,Vs_upper], [Vbulk_lower, Vp_lower, Vs_lower], rkkh, rkgh
            else:
                return [Vbulk_medium,Vp_medium,Vs_medium], [Vbulk_upper,Vp_upper,Vs_upper], [Vbulk_lower, Vp_lower, Vs_lower]

    def set_modal_composition(self, phase_list, fraction_list):
        
        """
        Prepares the input for HashinShtrikmanBounds method for the entered phase_list
        and fraction_list.
        
        Parameters:        
        - phase_list - list of phase ids to be included in the rock. 
            Example: ['fo', 'en'] -- Forsterite and Enstatite
        
        - fraction_list - list of fraction lists to be included in the rock
            Example: [0.8, 0.2] -- 0.8 fo and 0.2 en
            
        Returns:
        
        - phase_constant_list: 
            list containing the phase constats that can be inputted into HashinShtrikmanBounds method.
        - fraction_list: 
            np.ndarray list of fractions of entered minerals.
        """
        
        if type(fraction_list) != np.ndarray:
            fraction_list = np.array(fraction_list)
            
        all_sums = np.allclose(fraction_list.sum(axis = -1), 1.0, atol = 1e-6)
        
        if all_sums == False:
            raise ValueError('The values entered in fraction_list do not sum up to 1.')
        
        phase_constant_list = []
        
        for item in phase_list:
        
            phase_constant = self.get_phase_constants(phase = item)
            phase_constant_list.append(phase_constant)
                        
        return phase_constant_list, fraction_list
        
        