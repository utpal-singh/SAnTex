import pandas as pd
import numpy as np
from joblib import Parallel, delayed

def euler_to_quaternion(phi, theta, psi):
    """
    Convert Euler angles (in radians) to a quaternion.

    Parameters:
        phi (float): Euler angle phi in radians.
        theta (float): Euler angle theta in radians.
        psi (float): Euler angle psi in radians.

    Returns:
        np.ndarray: A quaternion representing the input Euler angles.    
    """
    c1 = np.cos(phi / 2)
    s1 = np.sin(phi / 2)
    c2 = np.cos(theta / 2)
    s2 = np.sin(theta / 2)
    c3 = np.cos(psi / 2)
    s3 = np.sin(psi / 2)

    q0 = c1 * c2 * c3 + s1 * s2 * s3
    q1 = s1 * c2 * c3 - c1 * s2 * s3
    q2 = c1 * s2 * c3 + s1 * c2 * s3
    q3 = c1 * c2 * s3 - s1 * s2 * c3

    return np.array([q0, q1, q2, q3])

def quaternion_to_angle(q1, q2):
    """
    Calculate the misorientation angle between two quaternions.

    Parameters:
        q1 (np.ndarray): First quaternion.
        q2 (np.ndarray): Second quaternion.

    Returns:
        float: The misorientation angle between the two quaternions in degrees.
    """
    dot_product = np.dot(q1, q2)
    angle_radians = 2 * np.arccos(np.abs(dot_product))
    return np.degrees(angle_radians)

def misorientation_angle(euler1, euler2):
    """
    Calculate the misorientation angle between two sets of Euler angles.

    Parameters:
        euler1 (list or np.ndarray): First set of Euler angles in radians.
        euler2 (list or np.ndarray): Second set of Euler angles in radians.

    Returns:
        float: The misorientation angle between the two sets of Euler angles in degrees.    
    """
    q1 = euler_to_quaternion(*np.radians(euler1))
    q2 = euler_to_quaternion(*np.radians(euler2))
    return quaternion_to_angle(q1, q2)

def assign_to_grains_parallel(df, threshold):
    """
        Assign grains to rows in a DataFrame based on misorientation angle threshold.

    Parameters:
        df (pd.DataFrame): Input DataFrame containing Euler angles for each row.
        threshold (float): Misorientation angle threshold in degrees.

    Returns:
        dict: A dictionary mapping row indices to assigned grain indices.
    """
    grains = []
    grain_indices = {}

    total_rows = len(df)
    processed_rows = 0

    def process_row(row):
        nonlocal grains, grain_indices, processed_rows
        euler_angles = row[['Euler1', 'Euler2', 'Euler3']].values
        assigned = False
        for grain_idx, grain_euler in enumerate(grains):
            misorientation = misorientation_angle(euler_angles, grain_euler)
            if misorientation <= threshold:
                grain_indices[row.name] = grain_idx
                grains[grain_idx] = np.mean([grain_euler, euler_angles], axis=0)
                assigned = True
                break
        if not assigned:
            grain_indices[row.name] = len(grains)
            grains.append(euler_angles)

        processed_rows += 1
        progress = (processed_rows / total_rows) * 100
        print(f"Progress: {progress:.2f}%\r", end='')

    Parallel(n_jobs=-1, prefer="threads")(delayed(process_row)(row) for _, row in df.iterrows())

    return grain_indices

