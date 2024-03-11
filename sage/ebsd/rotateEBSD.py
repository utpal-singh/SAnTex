import numpy as np
from scipy.spatial.transform import Rotation

def bunge_euler_rotation(phi1, Phi, phi2, angles):
    """
    Calculates the resultant Bunge Euler angles after applying a custom rotation.

    Args:
        phi1 (float): Initial Bunge Euler angle phi1 in degrees.
        Phi (float): Initial Bunge Euler angle Phi in degrees.
        phi2 (float): Initial Bunge Euler angle phi2 in degrees.
        angles (tuple): Custom rotation angles in degrees (alpha, beta, gamma).

    Returns:
        tuple: Resultant Bunge Euler angles (phi1, Phi, phi2) in degrees in the range [0, 360).
    """
    # Convert initial Bunge Euler angles to radians
    phi1_rad = np.deg2rad(phi1)
    Phi_rad = np.deg2rad(Phi)
    phi2_rad = np.deg2rad(phi2)

    # Create the initial rotation matrix from Bunge Euler angles
    initial_rot = Rotation.from_euler('zxz', [phi1_rad, Phi_rad, phi2_rad], degrees=False)

    # Create the custom rotation matrix
    alpha_rad, beta_rad, gamma_rad = [np.deg2rad(angle) for angle in angles]
    custom_rot = Rotation.from_euler('zxz', [alpha_rad, beta_rad, gamma_rad], degrees=False)

    # Calculate the resultant rotation matrix
    resultant_rot = initial_rot * custom_rot

    # Convert the resultant rotation matrix to Bunge Euler angles
    phi1_res, Phi_res, phi2_res = resultant_rot.as_euler('zxz', degrees=True)

    # Keep the angles in the range [0, 360)
    phi1_res = (phi1_res + 360) % 360
    Phi_res = (Phi_res + 360) % 360
    phi2_res = (phi2_res + 360) % 360

    return phi1_res, Phi_res, phi2_res

if __name__ == "__main__":
    # Example usage
    initial_phi1 = 119.468
    initial_Phi = 117.433
    initial_phi2 = 74.0688

    # Case 1: rotation('Euler', 180*degree, 180*degree, 0*degree, 'BUNGE,ZXZ')
    angles1 = (180, 180, 0)
    res_phi1, res_Phi, res_phi2 = bunge_euler_rotation(initial_phi1, initial_Phi, initial_phi2, angles1)
    print(f"Resultant Bunge Euler angles (degrees): phi1 = {res_phi1:.3f}, Phi = {res_Phi:.3f}, phi2 = {res_phi2:.3f}")

    # Case 2: rotation('Euler', 180*degree, 0*degree, 0*degree, 'BUNGE,ZXZ')
    angles2 = (180, 0, 0)
    res_phi1, res_Phi, res_phi2 = bunge_euler_rotation(initial_phi1, initial_Phi, initial_phi2, angles2)
    print(f"Resultant Bunge Euler angles (degrees): phi1 = {res_phi1:.3f}, Phi = {res_Phi:.3f}, phi2 = {res_phi2:.3f}")

    # Case 3: rotation('Euler', 180*degree, 180*degree, 180*degree, 'BUNGE,ZXZ')
    angles3 = (180, 180, 180)
    res_phi1, res_Phi, res_phi2 = bunge_euler_rotation(initial_phi1, initial_Phi, initial_phi2, angles3)
    print(f"Resultant Bunge Euler angles (degrees): phi1 = {res_phi1:.3f}, Phi = {res_Phi:.3f}, phi2 = {res_phi2:.3f}")