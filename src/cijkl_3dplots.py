import numpy as np
import math

import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

def christoffel_tensor(cijkl, n):
    tik = np.zeros((3, 3))

    for i in range(3):
        for k in range(3):
            for j in range(3):
                for l in range(3):
                    tik[i][k] += cijkl[i][j][k][l] * n[j] * n[l]
                    
    return tik

def wave_property(tik):
    eigenvalues, eigenvectors = np.linalg.eig(tik)
    indices = np.argsort(eigenvalues)[::-1]
    wave_moduli = []
    polarization_directions = []

    for i in indices:
        wave_moduli.append(eigenvalues[i])
        polarization_directions.append(eigenvectors[:, i] / np.linalg.norm(eigenvectors[:, i]))
    return tuple(wave_moduli), tuple(polarization_directions)

def phase_velocity(cijkl, rho):
    vp = []
    vs1 = []
    vs2 = []
    step = math.pi / 180

    for theta in np.arange(0, math.pi + step, step):
        for phi in np.arange(0, 2 * math.pi + step, step):
            n = np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])
            tik = christoffel_tensor(cijkl, n)

            wave_moduli, polarization_directions = wave_property(tik)

            vp.append(math.sqrt(wave_moduli[0] / rho))
            vs1.append(math.sqrt(wave_moduli[1] / rho))
            vs2.append(math.sqrt(wave_moduli[2] / rho))
    return tuple(vp), tuple(vs1), tuple(vs2)


from mpl_toolkits.mplot3d import Axes3D

def plot_phase_velocity(cijkl, rho):
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    step = math.pi / 180
    x = []
    y = []
    z = []
    c = []

    for theta in np.arange(0, math.pi + step, step):
        for phi in np.arange(0, 2 * math.pi + step, step):
            n = np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])
            tik = christoffel_tensor(cijkl, n)

            wave_moduli, polarization_directions = wave_property(tik)
            vp = math.sqrt(wave_moduli[0] / rho)
            vs1 = math.sqrt(wave_moduli[1] / rho)
            vs2 = math.sqrt(wave_moduli[2] / rho)

            x.append(n[0])
            y.append(n[1])
            z.append(n[2])

            c.append((vp + vs1 + vs2) / 3)

    ax.scatter(x, y, z, c=c)
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title('Phase velocity in 3d plots with velocity represented by colour')
    plt.show()

import numpy as np
import tensor_conversion
M = np.array([[198.96,   73.595,  68.185,   0.,      9.735,   0.   ],
                [ 73.595, 155.94,   62.23,    0.,      6.295,   0.   ],
                [ 68.185,  62.23,  225.99,    0.,     33.85,    0.   ],
                [  0.,      0.,      0.,     65.66,    0.,      6.415],
                [  9.735,   6.295,  33.85,    0.,     60.23,    0.   ],
                [  0.,      0.,     0.,      6.415,   0.,     65.18 ]])*10**9

cijkl = tensor_conversion.voigt_to_tensor(M)

rho = 3500

plot_phase_velocity(cijkl, rho)
