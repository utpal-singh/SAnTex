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
            #propagation direction vector n from theta and phi using spherical coordinates
            n = np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])
            # print(n)
            tik = christoffel_tensor(cijkl, n)
            wave_moduli, polarization_directions = wave_property(tik)
            vp.append(math.sqrt(wave_moduli[0] / rho))
            vs1.append(math.sqrt(wave_moduli[1] / rho))
            vs2.append(math.sqrt(wave_moduli[2] / rho))
    return tuple(vp), tuple(vs1), tuple(vs2)


def plot_vp_2d(cijkl, rho):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    step = math.pi / 180
    x = []
    y = []
    c = []

    for theta in np.arange(0, math.pi / 2 + step, step):
        for phi in np.arange(0, 2 * math.pi + step, step):
            n = np.array([math.sin(theta) * math.cos(phi), math.sin(theta) * math.sin(phi), math.cos(theta)])
            # print(n)
            tik = christoffel_tensor(cijkl, n)
            wave_moduli, polarization_directions = wave_property(tik)
            vp = math.sqrt(wave_moduli[0] / rho)
            # print(vp)
            x.append(n[0] / (1 + n[2]))
            y.append(n[1] / (1 + n[2]))
            c.append(vp)

    sc = ax.scatter(x, y, c=c, cmap='RdBu')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_title('vp in 2d stereographic projection plots with velocity represented by colour')
    cb = plt.colorbar(sc)
    cb.set_label('vp')
    cb.ax.set_yticklabels(['{:.1f}'.format(v) for v in cb.get_ticks()])
    plt.show()

if __name__=="__main__":
    from sage import Tensor
    tensor = Tensor()
    M = np.array([[198.96,   73.595,  68.185,   0.,      9.735,   0.   ],
                    [ 73.595, 155.94,   62.23,    0.,      6.295,   0.   ],
                    [ 68.185,  62.23,  225.99,    0.,     33.85,    0.   ],
                    [  0.,      0.,      0.,     65.66,    0.,      6.415],
                    [  9.735,   6.295,  33.85,    0.,     60.23,    0.   ],
                    [  0.,      0.,     0.,      6.415,   0.,     65.18 ]])*10**9

    cijkl = tensor.voigt_to_tensor(M)

    rho = 3500

    plot_vp_2d(cijkl, rho)
