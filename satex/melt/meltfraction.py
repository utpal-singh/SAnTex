def calculate_melt_fractions(pressures_input, temperatures_input, cutoff=1):
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

if __name__ == "__main__":
    pressure = [4, 5]  # GPa
    temperature = [2000, 3000]  # Kelvin
    cutoff_fraction = 0.05  # Cutoff value for melt fraction

    melt_fraction = calculate_melt_fractions(pressure, temperature, 0.05)
    print(melt_fraction)

    pressure = 4  # GPa
    temperature = 2000  # Kelvin
    cutoff_fraction = 0.05  # Cutoff value for melt fraction

    melt_fraction = calculate_melt_fractions(pressure, temperature, 0.05)
    print(melt_fraction)