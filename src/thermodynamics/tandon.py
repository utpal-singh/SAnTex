import numpy as np

def MS_tandon_and_weng(vp, vs, rho, vpi, vsi, rhoi, delta, c):
    # Conversion factor
    to_m_s = 1e3

    # Convert to m/s
    vp *= to_m_s
    vs *= to_m_s
    vpi *= to_m_s
    vsi *= to_m_s

    # Weighted average density
    rh = (1.0 - c) * rho + c * rhoi

    amu = vs ** 2 * rho
    amui = vsi ** 2 * rhoi
    alam = vp ** 2 * rho - 2.0 * amu
    alami = vpi ** 2 * rhoi - 2.0 * amui
    bmi = alami + amui * 2.0 / 3.0
    bmps = alam + amu

    # Young's modulus for matrix
    E0 = amu * (3.0 * alam + 2.0 * amu) / (alam + amu)
    # Poisson's ratio of the matrix.
    anu = alam / (2.0 * (alam + amu))

    # Some time-saving terms
    t1 = delta ** 2 - 1.0
    t2 = 1.0 - anu
    t3 = 1.0 - 2.0 * anu
    t4 = 3.0 * delta ** 2
    t5 = 1.0 - delta ** 2

    # D1, D2, and D3 from Tandon and Weng (1984)
    D1 = 1.0 + 2.0 * (amui - amu) / (alami - alam)
    D2 = (alam + 2.0 * amu) / (alami - alam)
    D3 = alam / (alami - alam)

    # g and g' terms
    if delta >= 1:
        acshdelta = np.log(delta + np.sqrt(t1))
        g = (delta * np.sqrt(t1) - acshdelta) * delta / (np.sqrt(t1) ** 3)
    else:
        # g' below
        g = (np.arccos(delta) - delta * np.sqrt(t5)) * delta / (np.sqrt(t5) ** 3)

    # Eshelby's Sijkl tensor
    s11 = (t3 + (t4 - 1.0) / t1 - (t3 + t4 / t1) * g) / (2.0 * t2)
    s22 = (t4 / (t1 * 2.0) + (t3 - 9.0 / (4.0 * t1)) * g) / (4.0 * t2)
    s33 = s22
    s23 = (delta ** 2 / (2.0 * t1) - (t3 + 3.0 / (4.0 * t1)) * g) / (4.0 * t2)
    s32 = s23
    s21 = (-2.0 * delta ** 2 / t1 + (t4 / t1 - t3) * g) / (4.0 * t2)
    s31 = s21
    s12 = (-1.0 * (t3 + 1.0 / t1) + (t3 + 3.0 / (2.0 * t1)) * g) / (2.0 * t2)
    s13 = s12
    s44 = (delta ** 2 / (2.0 * t1) + (t3 - 3.0 / (4.0 * t1)) * g) / (4.0 * t2)
    s66 = (t3 - (t1 + 2.0) / t1 - (t3 - 3.0 * (t1 + 2.0) / t1) * g / 2.0) / (4.0 * t2)
    s55 = s66

    # Tandon and Weng's B terms
    B1 = c * D1 + D2 + (1.0 - c) * (D1 * s11 + 2.0 * s21)
    B2 = c + D3 + (1.0 - c) * (D1 * s12 + s22 + s23)
    B3 = c + D3 + (1.0 - c) * (s11 + (1.0 + D1) * s21)
    B4 = c * D1 + D2 + (1.0 - c) * (s12 + D1 * s22 + s23)
    B5 = c + D3 + (1.0 - c) * (s12 + s22 + D1 * s23)

    # Tandon and Weng's A terms
    A1 = D1 * (B4 + B5) - 2.0 * B2
    A2 = (1.0 + D1) * B2 - (B4 + B5)
    A3 = B1 - D1 * B3
    A4 = (1.0 + D1) * B1 - 2.0 * B3
    A5 = (1.0 - D1) / (B4 - B5)
    A = 2.0 * B2 * B3 - B1 * (B4 + B5)

    # Tandon and Weng (1984) equations (25), (28), (31), (32)
    E11 = E0 / (1.0 + c * (A1 + 2.0 * anu * A2) / A)
    E22 = E0 / (1.0 + c * (-2.0 * anu * A3 + (1.0 - anu) * A4 + (1.0 + anu) * A5 * A) / (2.0 * A))
    amu12 = amu * (1.0 + c / (amu / (amui - amu) + 2.0 * (1.0 - c) * s66))
    amu23 = amu * (1.0 + c / (amu / (amui - amu) + 2.0 * (1.0 - c) * s44))

    # Sayers equation (36)
    anu31 = anu - c * (anu * (A1 + 2.0 * anu * A2) + (A3 - anu * A4)) / (A + c * (A1 + 2.0 * anu * A2))

    # T&W equation (36)
    # aK12 term; bmps=plane strain bulk modulus
    anum = (1.0 + anu) * (1.0 - 2.0 * anu)
    denom = 1.0 - anu * (1.0 + 2.0 * anu31) + c * (2.0 * (anu31 - anu) * A3 + (1.0 - anu * (1.0 + 2.0 * anu31)) * A4) / A
    aK23 = bmps * anum / denom
    anu12tst = E11 / E22 - (1.0 / amu23 + 1.0 / aK23) * E11 / 4.0

    # Cij - Sayers' (1992) equations (24)-(29)
    # Conversion
    CC = np.zeros((6, 6))
    CC[1, 1] = amu23 + aK23
    CC[2, 2] = CC[1, 1]
    CC[3, 3] = CC[1, 1]
    CC[4, 4] = (CC[2, 2] - CC[2, 3]) / 2.0

    # Fill out matrix by symmetry
    for i in range(6):
        for j in range(i, 6):
            CC[j, i] = CC[i, j]

    # Convert to GPa
    CC /= 1e9

    return CC

# Example usage:
vp = 3000.0  # P-wave velocity in m/s
vs = 1500.0  # S-wave velocity in m/s
rho = 2000.0  # Density in kg/m^3
vpi = 4000.0  # P-wave velocity of inclusion in m/s
vsi = 2000.0  # S-wave velocity of inclusion in m/s
rhoi = 2500.0  # Density of inclusion in kg/m^3
delta = 1.5  # Aspect ratio of inclusion
c = 0.2  # Volume fraction of inclusion

result = MS_tandon_and_weng(vp, vs, rho, vpi, vsi, rhoi, delta, c)
print(result)
