from satex import Anisotropy
import numpy as np

if __name__ == "__main__":
    # Example stiffness matrix and density
    stiffness_matrix = np.array([[198.96, 73.595, 68.185, 0., 9.735, 0.],
                                [73.595, 155.94, 62.23, 0., 6.295, 0.],
                                [68.185, 62.23, 225.99, 0., 33.85, 0.],
                                [0., 0., 0., 65.66, 0., 6.415],
                                [9.735, 6.295, 33.85, 0., 60.23, 0.],
                                [0., 0., 0., 6.415, 0., 65.18]]) * 10**9

    density = 3500

    # Create an instance of the Anisotropy class
    anisotropy_instance = Anisotropy(stiffness_matrix, density)
    anisotropy_instance.plot()