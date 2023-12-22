import math

class ElasticConstant:
    def __init__(self, E=9, n=0.01, G=421, K=3, l=2):
        self.E = E  # Young's modulus
        self.n = n  # Poisson's ratio
        self.G = G  # Shear modulus
        self.K = K  # Bulk modulus
        self.l = l  # Lamé parameter

        # Check value to satisfy thermodynamic 2nd law limit: info on https://www.efunda.com/formulae/solid_mechanics/mat_mechanics/calc_elastic_constants.cfm#calc
        if E <= 0 or G <= 0 or K <= 0:
            raise ValueError("E, G and K should be greater than 0")
        if n < -1 or n > 0.5:
            raise ValueError("n should be between -1 and 0.5 inclusive")

        # Calculating R for E, l input case
        self.R = math.sqrt(E**2 +9*l**2 +2*E*l)

    def en(self):
        # For E, n input case
        G_en = self.E / (2 * (1 + self.n))
        K_en = self.E / (3 * (1 - 2*self.n))
        l_en = 2 * G_en * self.n / (1 - 2*self.n)

        # Check the output to satisfy the limit
        if G_en <= 0 or K_en <= 0 or l_en <= 0:
            raise ValueError("G, K and l must be greater than 0")

        return G_en, K_en, l_en

    def eg(self):
        # For E, G input case
        n_eg = self.E / (2*self.G) -1 
        K_eg =self.E * self.G / (9*self.G-3*self.E)
        l_eg =(3*K_eg*(3*K_eg-self.E))/(9*K_eg-self.E)

        # Check the output to satisfy the limit
        if n_eg < -1 or n_eg > 0.5 or K_eg <= 0 or l_eg <= 0:
            raise ValueError("n must be between -1 and 0.5 inclusive, K and l must be greater than 0")

        return n_eg, K_eg, l_eg

    def ek(self):
        # For E, K input case
        n_ek = (3*self.K - self.E) / (6*self.K)
        G_ek = (3*self.K*self.E) / (9*self.K - self.E)
        l_ek = 3*self.K * (3*self.K - self.E) / (9*self.K - self.E)

        # Check the output to satisfy the limit
        if n_ek < -1 or n_ek > 0.5 or G_ek <= 0 or l_ek <= 0:
            raise ValueError("n must be between -1 and 0.5 inclusive, G and l must be greater than 0")

        return n_ek, G_ek, l_ek

    def el(self):
        # For E, l input case
        n_el = 2*self.l / (self.E + self.l + self.R)
        G_el = (self.E - 3*self.l + self.R) / 4
        K_el = (self.E + 3*self.l + self.R) / 6

        # Check the output to satisfy the limit
        if n_el < -1 or n_el > 0.5 or G_el <= 0 or K_el <= 0:
            raise ValueError("n must be between -1 and 0.5 inclusive, G and K must be greater than 0")

        return n_el, G_el, K_el
    
    def ng(self):
        e_ng = 2*self.G*(1 + self.n)
        k_ng = 2*self.G*(1+self.n)/(3*(1-2*self.n))
        l_ng = 2*self.G*self.n/(1-2*self.n)

        return e_ng, k_ng, l_ng

    def nk(self):
        e_nk = 3*self.K*(1 - 2*self.n)
        g_nk = 3*self.K*(1-2*self.n)/(2*(1+ self.n))
        l_nk = 3*self.K*self.n/(1+self.n)
        return e_nk, g_nk, l_nk

    def nl(self):
        e_nl = l(1+n)(1-2n)/n
        g_nl = l(1-2n)/2(1+n)
        k_nl = l(1+n)/3n
        return e_nl, g_nl, k_nl

    def gk(self):
        # For G, K input case
        E_gk = (9*self.K*self.G) / (3*self.K + self.G)
        n_gk = (3*self.K - 2*self.G) / (2*(3*self.K + self.G))
        l_gk = (3*self.K - 2*self.G) / 3

        # Check the output to satisfy the limit
        if E_gk <= 0 or n_gk < -1 or n_gk > 0.5 or l_gk <= 0:
            raise ValueError("E and l must be greater than 0, n must be between -1 and 0.5 inclusive")

        return E_gk, n_gk, l_gk

    def gl(self):
        # For G, l input case
        E_gl = self.R - 2*self.G - self.l
        n_gl = (self.R - 2*self.G - self.l) / (4*self.G)
        K_gl = (self.R + 2*self.G + self.l) / 3

        # Check the output to satisfy the limit
        if E_gl <= 0 or n_gl < -1 or n_gl > 0.5 or K_gl <= 0:
            raise ValueError("E and K must be greater than 0, n must be between -1 and 0.5 inclusive")

        return E_gl, n_gl, K_gl

    def kl(self):
        # For K, l input case
        E_kl = (3*self.K - self.l) / (1 - self.n)
        n_kl = (3*self.K - self.l - self.R) / (6*self.K - 2*self.l)
        G_kl = (3*self.K - self.l - self.R) / 2

        # Check the output to satisfy the limit
        if E_kl <= 0 or n_kl < -1 or n_kl > 0.5 or G_kl <= 0:
            raise ValueError("E and G must be greater than 0, n must be between -1 and 0.5 inclusive")

        return E_kl, n_kl, G_kl
    

if __name__ == "__main__":
    elastic = ElasticConstant(n = 0.08, K = 37.6)
    a, b, c = elastic.nk()
    print("e: ", a)
    print("g: ", b)
    print("l: ", c)
