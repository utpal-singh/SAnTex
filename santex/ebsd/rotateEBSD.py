import numpy as np
from scipy.spatial.transform import Rotation
from joblib import Parallel, delayed
from tqdm import tqdm

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

def apply_custom_rotation_to_dataframe(euler_df, angles):
    """
    Applies a custom rotation to each row of a DataFrame containing Bunge Euler angles.

    Args:
        euler_df (pd.DataFrame): DataFrame containing Bunge Euler angles.
        angles (tuple): Custom rotation angles in degrees (alpha, beta, gamma).

    Returns:
        pd.DataFrame: DataFrame with updated Bunge Euler angles after rotation.
    """
    def bunge_euler_rotation_parallel(row):
        """
        Calculates the resultant Bunge Euler angles after applying a custom rotation in parallel.

        Args:
            row (pd.Series): Row of DataFrame containing Euler angles.

        Returns:
            tuple: Resultant Bunge Euler angles (phi1, Phi, phi2) in degrees in the range [0, 360).
        """
        # Extract Euler angles from the row
        phi1, Phi, phi2 = row['Euler1'], row['Euler2'], row['Euler3']
        
        # Call the original rotation function
        return bunge_euler_rotation(phi1, Phi, phi2, angles)

    # Calculate the number of iterations
    num_iterations = len(euler_df)

    # Parallelize the process and show progress in percent
    results = Parallel(n_jobs=-1)(
        delayed(bunge_euler_rotation_parallel)(row)
        for _, row in tqdm(euler_df.iterrows(), total=num_iterations)
    )

    # Update the DataFrame with the new Euler angles
    euler_df[['Euler1', 'Euler2', 'Euler3']] = results

    return euler_df

def apply_custom_rotation_to_dataframe_noxy(df, angles):
    """
    Applies a custom rotation to each row of a DataFrame containing Bunge Euler angles, X, and Y coordinates.

    Args:
        df (pd.DataFrame): DataFrame containing Bunge Euler angles, X, and Y coordinates.
        angles (tuple): Custom rotation angles in degrees (alpha, beta, gamma).

    Returns:
        pd.DataFrame: DataFrame with updated Bunge Euler angles, X, and Y coordinates after rotation.
    """
    def rotation_parallel(row):
        """
        Calculates the resultant Bunge Euler angles, X, and Y coordinates after applying a custom rotation in parallel.

        Args:
            row (pd.Series): Row of DataFrame containing Euler angles, X, and Y coordinates.

        Returns:
            tuple: Resultant Bunge Euler angles (phi1, Phi, phi2) in degrees in the range [0, 360),
                   rotated X and Y coordinates.
        """
        # Extract Euler angles, X, and Y coordinates from the row
        phi1, Phi, phi2 = row['Euler1'], row['Euler2'], row['Euler3']
        x, y = row['X'], row['Y']
        
        # Apply rotation to Euler angles
        phi1_res, Phi_res, phi2_res = bunge_euler_rotation(phi1, Phi, phi2, angles)
        
        # Apply rotation to X and Y coordinates
        rotated_coords = Rotation.from_euler('zxz', [np.deg2rad(angle) for angle in angles], degrees=False).apply([x, y, 0])
        
        return phi1_res, Phi_res, phi2_res, rotated_coords[0], rotated_coords[1]

    # Calculate the number of iterations
    num_iterations = len(df)

    # Parallelize the process and show progress in percent
    results = Parallel(n_jobs=-1)(
        delayed(rotation_parallel)(row)
        for _, row in tqdm(df.iterrows(), total=num_iterations)
    )

    # Update the DataFrame with the new Euler angles, X, and Y coordinates
    df[['Euler1', 'Euler2', 'Euler3', 'X', 'Y']] = results

    return df
