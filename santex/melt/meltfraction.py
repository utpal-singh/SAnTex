import numpy as np
def calculate_melt_fractions(pressures_input, temperatures_input, cutoff=1):
    """
    Calculates melt fractions based on pressure and temperature conditions.

    Parameters:
    pressures_input (float, list, or ndarray): Pressure or list/array of pressures in GPa.
    temperatures_input (float, list, or ndarray): Temperature or list/array of temperatures in Kelvin.
    cutoff (float, optional): Cutoff value for melt fraction calculation. Default is 1.

    Returns:
    float or ndarray: Melt fraction(s) calculated based on the input conditions.

    Notes:
    - The function computes melt fractions based on pressure and temperature using thermodynamic equations.
    - If pressures_input and temperatures_input are arrays of the same shape, returns an array of melt fractions.
    - If pressures_input and temperatures_input are scalars, returns a single melt fraction.
    - Values outside the cutoff range are capped at the cutoff value.
    """
    pressures = np.atleast_1d(pressures_input)
    temperatures = np.atleast_1d(temperatures_input)
    
    if pressures.shape != temperatures.shape:
        raise ValueError("The pressures and temperatures must have compatible shapes.")
    
    T_solidus = np.where(pressures <= 15,
                         1358.7 + 132.899 * pressures - 5.104 * pressures**2,
                         1510.76 + 46.27 * pressures - 0.8036 * pressures**2)
    T_liquidus = np.where(pressures <= 15,
                          2053.0 + 45.0 * pressures - 2.0 * pressures**2,
                          1470.3025 + 55.53 * pressures - 0.9084 * pressures**2)

    T_ss = (temperatures - (T_solidus + T_liquidus) / 2) / (T_liquidus - T_solidus)
    F = np.where((-0.5 <= T_ss) & (T_ss <= 0.5),
                 0.5 + T_ss + (T_ss**2 - 0.25) * (0.4256 + 2.988 * T_ss),
                 np.where(T_ss < -0.5, 0, 1))
    
    F = np.where((F > cutoff) , cutoff, F)
    
    return F if isinstance(pressures_input, np.ndarray) or isinstance(pressures_input, list) else F.item()
