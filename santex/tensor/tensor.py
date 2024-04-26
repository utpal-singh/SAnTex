import numpy as np
class Tensor:
    def __init__(self, data=None):

        """
        Initialize a Tensor object.

        Parameters:
        - data (array): The data to initialize the tensor with. If None, an empty tensor is created.

        Returns:
        - None
        """
        self.data = np.array(data) if data is not None else None

    def tensor_to_voigt(self, tensor):

        """
        Convert a tensor to its Voigt notation.

        Parameters:
        - tensor (array): The input tensor.

        Returns:
        - array-like: The tensor converted to Voigt notation.
        """
        voigt = [[0 for i in range(6)] for j in range(6)]
        mapping = {(1, 1): 1, (2, 2): 2, (3, 3): 3, (2, 3): 4, (3, 2): 4, (1, 3): 5, (3, 1): 5, (1, 2): 6, (2, 1): 6}

        for i in range(1, 4):
            for j in range(1, 4):
                for k in range(1, 4):
                    for l in range(1, 4):

                        I = mapping[(i, j)]
                        J = mapping[(k, l)]
                        voigt[I-1][J-1] = tensor[i-1][j-1][k-1][l-1]

        return np.array(voigt)


    def voigt_to_tensor(self, voigt):
        """
        Convert Voigt notation back to a tensor.

        Parameters:
        - voigt (array): The input Voigt notation.

        Returns:
        - array: The Voigt notation converted back to a tensor.
        """
        tensor = [[[[0 for i in range(3)] for j in range(3)] for k in range(3)] for l in range(3)]
        mapping = {1: (1, 1), 2: (2, 2), 3: (3, 3), 4: (2, 3), 5: (1, 3), 6: (1, 2)}
        
        for I in range(1, 7):
            for J in range(1, 7):
                i, j = mapping[I]
                k, l = mapping[J]
                tensor[i-1][j-1][k-1][l-1] = voigt[I-1][J-1]
                tensor[j-1][i-1][k-1][l-1] = voigt[I-1][J-1]
                tensor[i-1][j-1][l-1][k-1] = voigt[I-1][J-1]
                tensor[j-1][i-1][l-1][k-1] = voigt[I-1][J-1]
        return np.array(tensor)
    
    def euler_to_rotation(self, phi1, phi, phi2):

        """
        Convert Euler angles to a rotation matrix.

        Parameters:
        - phi1 (float): First Euler angle in degrees.
        - phi (float): Second Euler angle in degrees.
        - phi2 (float): Third Euler angle in degrees.

        Returns:
        - array: The rotation matrix corresponding to the Euler angles.
        """
        phi1 = np.radians(phi1)
        phi = np.radians(phi)
        phi2 = np.radians(phi2)
        Rz1 = np.array([[np.cos(phi1), -np.sin(phi1), 0], [np.sin(phi1), np.cos(phi1), 0], [0, 0, 1]])
        Rx = np.array([[1, 0, 0], [0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi)]])
        Rz2 = np.array([[np.cos(phi2), -np.sin(phi2), 0], [np.sin(phi2), np.cos(phi2), 0], [0, 0, 1]])
        R = Rz1 @ Rx @ Rz2 

        return R

    def rotate_tensor(self, tensor, phi1, phi, phi2):
        """
        Rotate a tensor using Euler angles.

        Parameters:
        - tensor (array-like): The input tensor to be rotated.
        - phi1 (float): First Euler angle in degrees.
        - phi (float): Second Euler angle in degrees.
        - phi2 (float): Third Euler angle in degrees.

        Returns:
        - array: The rotated tensor.
        """
        R = self.euler_to_rotation(phi1, phi, phi2)
        tensor = np.array(tensor)
        R_tensor = np.einsum('im,jn,ko,lp,mnop->ijkl', R, R, R, R, tensor)

        return R_tensor

