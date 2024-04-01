import numpy as np

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
    
def solidus(P):
    if (P <= 15.0e9):
        return 1358.7 + 132.899*P/(1.0e9)- 5.104*(P/(1.0e9))**2
    else:
        return 1510.76 + 46.27*P/(1.0e9)- 0.8036*(P/(1.0e9))**2
    
# Solidus array
P_range = np.arange(0.0,9e9+1,1e6)
Solidus = []
for P in P_range:
    TSol = solidus(P)
    Solidus.append(TSol)

# Liquidus function
def liquidus(P):
    if (P <= 15.0e9):
        return 2053. + 45*P/(1.0e9)- 2*(P/(1.0e9))**2
    else:
        return 1470.3025 + 55.53*P/(1.0e9)- 0.9084*(P/(1.0e9))**2
    
# Liquidus array
P_range = np.arange(0.0,9e9+1,1e6)
Liquidus = []
for P in P_range:
    TLiq = liquidus(P)
    Liquidus.append(TLiq)

# Super solidus function
g=9.81
rho=(((2750.*40.)+(3310.*110.))/150.)
mid=(solidus(rho*g*100000)+liquidus(rho*g*100000))/2
delta=liquidus(rho*g*100000)-solidus(rho*g*100000)

def supersolidus(z):
    if (geotherm(z) > liquidus((rho*g*z)/1.0e9)): #P is expressed in Giga Pascal
        return 0.5
    elif (geotherm(z) > solidus((rho*g*z)/1.0e9)):
        return (geotherm(z)-mid)/delta
    else:
        return -0.5
    
# Super solidus array
P_range = np.arange(0.0,9e9+1,1e6)
SuperSolidus40l = []
SuperSolidus20l = []
SuperSolidus0 = []
SuperSolidus20r = []
SuperSolidus40r = []

for P in P_range:
    Tss40l = -0.4*(liquidus(P)-solidus(P))+(liquidus(P)+solidus(P))/2
    Tss20l = -0.2*(liquidus(P)-solidus(P))+(liquidus(P)+solidus(P))/2
    Tss0 = (liquidus(P)+solidus(P))/2
    Tss20r = 0.2*(liquidus(P)-solidus(P))+(liquidus(P)+solidus(P))/2
    Tss40r = 0.4*(liquidus(P)-solidus(P))+(liquidus(P)+solidus(P))/2
    SuperSolidus40l.append(Tss40l)    
    SuperSolidus20l.append(Tss20l)
    SuperSolidus0.append(Tss0)
    SuperSolidus20r.append(Tss20r)
    SuperSolidus40r.append(Tss40r)

# Melt fraction function
def meltfraction(z):
    return 0.5+supersolidus(z)+(supersolidus(z)**2-0.25)*(0.4256+2.988*supersolidus(z))

from sympy import symbols, solve

T_ss = symbols('T_ss')

F_T_ss_value = 0.8

equation = 0.5 + T_ss + (T_ss**2 - 0.25) * (0.4256 + 2.988 * T_ss) - F_T_ss_value

solutions = solve(equation, T_ss)

real_solutions = [sol.evalf() for sol in solutions if sol.is_real]
Tss = real_solutions
print(real_solutions)

def meltfraction(z):
    return 0.5+supersolidus(z)+(supersolidus(z)**2-0.25)*(0.4256+2.988*supersolidus(z))

print(supersolidus(400000))
print(liquidus(4e9))
print(geotherm(4e9))
print(meltfraction(270000))
print(0.5+(-0.69)+((-0.69)**2-0.25)*(0.4256+2.988*(-0.69)))
