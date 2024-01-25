import numpy as np

def euler_to_rotation(phi1, phi, phi2):
    phi1 = np.radians(phi1)
    phi = np.radians(phi)
    phi2 = np.radians(phi2)
    Rz1 = np.array([[np.cos(phi1), -np.sin(phi1), 0], [np.sin(phi1), np.cos(phi1), 0], [0, 0, 1]])
    Rx = np.array([[1, 0, 0], [0, np.cos(phi), -np.sin(phi)], [0, np.sin(phi), np.cos(phi)]])
    Rz2 = np.array([[np.cos(phi2), -np.sin(phi2), 0], [np.sin(phi2), np.cos(phi2), 0], [0, 0, 1]])
    R = Rz1 @ Rx @ Rz2 

    return R

def rotate_tensor(tensor, phi1, phi, phi2):
    R = euler_to_rotation(phi1, phi, phi2)
    tensor = np.array(tensor)
    R_tensor = np.einsum('im,jn,ko,lp,mnop->ijkl', R, R, R, R, tensor)

    return R_tensor
