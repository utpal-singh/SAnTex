Qm = 0.02  # Mantle heatflow in W per meter square
K = 2.25  # Conductivity in W/m/K

zc = 40e3 # Thickness of the crust in meters
zlab = 150e3  # Thickness of thermal lithosphere in meters
 
To = 293. # Surface temperature in Kelvin.
Tlab = 1603. # Temperature at LAB in Kelvin.

H = 0.518e-6 # Radiogenic heat production in the crust in Watt per cubic meter

Adiabat = 0.0003 # Adiabatic gradient: 0.3K per kilometer

# Using the geothem equation for the crust, the temperature at the Moho is given by:
Tzc = -(H/(2*K))*zc**2 + ((Qm/K) + (H/K)*zc)*zc + To

def geotherm(y):
    if (y  <= 0):
        return To
    elif (0 < y <= zc):
        return -(H/(2*K))*y**2+((Qm/K) + (H/K)*zc)*y + To
    elif (zc < y <= zlab):
        return Tzc + (Tlab - Tzc)*(y-zc)/(zlab-zc)
    else:
        return Tlab + 0.0003*(y - zlab)
    