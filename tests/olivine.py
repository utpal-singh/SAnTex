
import numpy as np
from scipy.optimize import fsolve

c = np.array([[323.70, 66.40, 71.60, .000, .000, .000], [66.40, 197.60, 
     75.60, .000, .000, .000], [71.60, 75.60, 
     235.10, .000, .000, .000], [.000, .000, .000, 
     64.62, .000, .000], [.000, .000, .000, .000, 
     78.05, .000], [.000, .000, .000, .000, .000, 79.04]]) * 10**9
phi=np.pi/4
th = np.pi/3
density=3.5
V= 10

l = np.cos(phi) * np.cos(np.pi/2 - th)
m = np.cos(np.pi/2 - th) * np.cos(np.pi/2 - phi)
n = np.cos(th)
A = c[0][0] * l**2 + c[5][5] * m**2 + c[4][4] * n**2 + 2 * c[0][5] * l * m + \
   2 * c[4][5] * m * n + 2 * c[0][4] * n * l
B = c[5][5] * l**2 + c[1][1] * m**2 + c[3][3] * n**2 + 2 * c[1][5] * l * m + \
   2 * c[1][3] * m * n + 2 * c[3][5] * n * l
EE = c[4][4] * l**2 + c[3][3] * m**2 + c[2][2] * n**2 + \
   2 * c[3][4] * l * m + 2 * c[2][3] * m * n + 2 * c[2][4] * n * l
F = c[4][5] * l**2 + c[1][3] * m**2 + \
   c[2][3] * n**2 + (c[1][4] + c[3][5]) * l * m + (c[1][2] + \
      c[3][3]) * m * n + (c[2][5] + c[3][4]) * n * l
G = c[0][4] * l**2 + c[3][5] * m**2 + \
   c[2][4] * n**2 + (c[0][3] + c[4][5]) * l * m + (c[2][5] + \
      c[3][4]) * m * n + (c[0][2] + c[4][4]) * n * l
H = c[0][5] * l**2 + c[1][5] * m**2 + \
   c[3][4] * n**2 + (c[0][1] + c[5][5]) * l * m + (c[1][4] + \
      c[3][5]) * m * n + (c[0][3] + c[4][5]) * n * l
velocity = np.array([[A - density*V**2, H, G], [H, B - density*V**2, F], [G, F,
    EE - density*V**2]])

def det_velocity(V):
    return np.linalg.det(velocity)

Ro = fsolve(det_velocity,V)
pwave = lambda th, phi: np.real(Ro[-1])
s1wave = lambda th, phi: np.real(Ro[-3])
s2wave = lambda th, phi: np.real(Ro[-6])
splitting = lambda th, phi: s1wave(th,phi) - s1wave(th,phi)