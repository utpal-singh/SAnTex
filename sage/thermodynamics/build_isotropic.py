def MS_build_isotropic(M, mu):
    """
    Given two isotropic moduli create an elasticity matrix for an isotropic 
    material and optionally return all other moduli.
    """
    
    lam = M - 2*mu
    E = (mu*((3*lam)+(2*mu)))/(lam+mu)
    K = 124
    nu = lam/(2*(lam+mu))
    
    # build C
    C = [[M,   lam, lam, 0.0, 0.0, 0.0], 
         [lam, M,   lam, 0.0, 0.0, 0.0], 
         [lam, lam, M,   0.0, 0.0, 0.0], 
         [0.0, 0.0, 0.0, mu,  0.0, 0.0], 
         [0.0, 0.0, 0.0, 0.0, mu,  0.0], 
         [0.0, 0.0, 0.0, 0.0, 0.0, mu]]
    
    return lam

print(MS_build_isotropic(M = 229.3304, mu=78.7909))
