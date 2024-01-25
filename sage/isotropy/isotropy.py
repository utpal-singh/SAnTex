import numpy as np
from scipy.optimize import root_scalar

def calc_seismic_prop(tt, pp, minprops, return_vp_vs_vbulk=False, return_aktout=False):
    """
    Calculate velocities at P,T for a mineral given parameters.
    """
    # Special handling if this is quartz
    if minprops['name'] in ['aqz', 'bqz', 'qz']:
        rho, aks, amu, akt = ah16_min_qz(tt, pp, minprops)
        if return_vp_vs_vbulk:
            vbulk, vp, vs = calculate_velocities(rho, aks, amu)
            if return_aktout:
                return rho, aks, amu, vp, vs, vbulk, akt
            return rho, aks, amu, vp, vs, vbulk
        return rho, aks, amu

    # Some numbers
    tref = 25.0
    tabs = 273.0
    tol = 1.e-8

    # input
    rho0 = minprops['rho0']
    ao = minprops['ao']
    akt0 = minprops['akt0']
    dkdp = minprops['dkdp']
    amu0 = minprops['amu0']
    dmudp = minprops['dmudp']
    gam = minprops['gam']
    grun = minprops['grun']
    delt = minprops['delt']

    # first, integrate with temperature at P=0
    alpht = ao * (1. - 10. / np.sqrt(tabs + tt))
    phi = ao * (tt - tref) - 20. * ao * (np.sqrt(tt + tabs) - np.sqrt(tref + tabs))
    rhot = rho0 * np.exp(-phi)

    aktt = akt0 * np.exp(-delt * phi)  # KtT in Excel macro
    amut = amu0 * np.exp(-gam * phi)   # GT in Excel macro

    # second, integrate with pressure at T fixed
    fguess = pp / (3. * aktt)
    fguess2 = 0.
    press = pp
    akk = aktt
    akkprime = dkdp

    res = root_scalar(pfunc, args=(press, akk, akkprime), bracket=[fguess, fguess2], xtol=tol)
    if not res.converged:
        raise ValueError("Root finding did not converge")
    
    ff = res.root
    dvol = (1. + 2. * ff) ** 1.5
    rho = rhot * dvol
    ffac = (1. + 2. * ff) ** 2.5

    akt = aktt * ffac * (1. - (5. - 3. * dkdp) * ff)
    amu = amut * ffac * (1. - ff * (5. - 3. * dmudp * aktt / amut))

    # convert to adiabatic bulk modulus
    alphtp = alpht * dvol ** (-delt)
    aks = akt * (1. + alphtp * grun * (tt + tabs))

    if return_vp_vs_vbulk:
        vbulk, vp, vs = calculate_velocities(rho, aks, amu)
        if return_aktout:
            return rho, aks, amu, vp, vs, vbulk, akt
        return rho, aks, amu, vp, vs, vbulk

    return rho, aks, amu

def calculate_velocities(rho, aks, amu):
    vbulk = 0.001 * np.sqrt(aks / rho)
    vp = 0.001 * np.sqrt((aks + 1.3333 * amu) / rho)
    vs = 0.001 * np.sqrt(amu / rho)
    return vbulk, vp, vs

def pfunc(f, press, akk, akkprime):
    """
    Pressure vs strain function used for minimization (root finding).
    """
    zeta = 0.75 * (4. - akkprime)
    f1 = (1. + 2. * f) ** 2.5
    f2 = (1. - 2. * zeta * f)
    pdif = 3. * akk * f * f1 * f2 - press
    return pdif

minprops = {
    'name': 'coe',
    'rho0': 2911,
    'ao': 1.8e-5,
    'akt0': 9.74e10,
    'dkdp': 4.3,
    'amu0': 6.16e10,
    'dmudp': 1.0541,
    'gam': 4.66,
    'grun': 0.36,
    'delt': 4.66,
}

rho, aks, amu, vp, vs, vbulk, akt = calc_seismic_prop(2000, 2e9, minprops, return_vp_vs_vbulk=True, return_aktout=True)
print(vp, vs)