import numpy as np

class Melt:
    def __init__(self, E, n):
        self.E = E
        self.n = n

    def get_elasticity_matrix(self):
        const_o = self.E/((1+self.n)*(1-2*self.n))
        elasticity_matrix_melt = np.zeros((6, 6))

        elasticity_matrix_melt[0,0] = 1- self.n
        elasticity_matrix_melt[1, 1] = 1- self.n
        elasticity_matrix_melt[2, 2] = 1- self.n
        elasticity_matrix_melt[3, 3] = 0.5 * (1 - 2*self.n)
        elasticity_matrix_melt[4, 4] = 0.5 * (1 - 2*self.n)
        elasticity_matrix_melt[5, 5] = 0.5 * (1 - 2*self.n)

        # print(const_o*(1- self.n), const_o*0.5*(1-2*self.n), const_o*self.n, 1/self.E)


if __name__=="__main__":
    from elastic_constants import ElasticConstant

    y = ElasticConstant(K = 16.1, G=0.01)
    e, n, l = y.gk()
    print(e)
    print(n)
    x = Melt(E = e, n = n)
    x.get_elasticity_matrix()
