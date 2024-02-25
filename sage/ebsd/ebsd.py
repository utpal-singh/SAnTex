# sage/ebsd/ebsd.py

from .ctf_parser import Ctf  # Importing Ctf class from ctf_parser module
import matplotlib.pyplot as plt  # Importing matplotlib for plotting functionalities

import numpy as np  # Importing numpy for numerical operations
import pandas as pd  # Importing pandas for data manipulation
from sklearn.neighbors import KNeighborsClassifier  # Importing KNeighborsClassifier for classification
from joblib import Parallel, delayed  # Importing Parallel and delayed for parallel execution
from ..tensor import Tensor
from scipy.spatial.transform import Rotation


class EBSD:
    """
    Class for handling Electron Backscatter Diffraction (EBSD) data.
    """

    def __init__(self, filename) -> None:
        """
        Initializes the EBSD object.

        Parameters:
            filename (str): The filename of the EBSD data file.
        """
        self._filename = filename
        self.ctf = Ctf(self._filename)  # Creating an instance of the Ctf class
    
    def get_euler_angles(self, phase, data = None):
        """
        Retrieves Euler angles for a specified phase.

        Parameters:
            phase (int): Integer representing the phase number.
            data (pandas.DataFrame, optional): DataFrame containing EBSD data. If not provided, uses stored data.


        Returns:
            pandas.DataFrame: DataFrame containing Euler angles (Euler1, Euler2, Euler3).
        """

        if data is None:
            data, header_data = self.ctf.get_data()
        phase_df = data[data['Phase'] == phase]
        euler_angles = phase_df[['Euler1', 'Euler2', 'Euler3']]
        euler_angles = euler_angles.reset_index(drop=True)
        return euler_angles
    
    def plot(self, data=None):
        """
        Plots the EBSD map with colors based on phase.

        Parameters:
            data (pandas.DataFrame, optional): DataFrame containing EBSD data. If not provided, uses stored data.

        Returns:
            None
        """
        if data is None:
            data, _ = self.ctf.get_data()
        df_new = data[data['Phase'].isin([0, 1, 2])]
        plt.scatter(df_new['X'], df_new['Y'], c=df_new['Phase'], cmap='viridis')
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.title('EBSD Map with Color Based on Phase')
        plt.show()

    def phases_names(self):
        """
        Retrieves phase names.

        Returns:
            pandas.DataFrame: DataFrame containing phase names.
        """
        phases_info =  self.ctf.phases_info()[6]
        column_keep = ["phase"]
        phases_ = phases_info[column_keep].copy()
        return phases_
    
    def phases_info(self):
        """
        Retrieves phase information.

        Returns:
            pandas.DataFrame: DataFrame containing phase information.
        """
        phases_info_ =  self.ctf.phases_info()[6]
        return phases_info_
    
    def get_euler_angles_by_phase(self, phase):
        """
        Retrieves Euler angles by phase names.
        Not yet implemented.

        Parameters:
            phase (str): Phase name.

        Returns:
            None
        """
        pass

    def get_phases_numbers(self):
        """
        Retrieves phase numbers.

        Returns:
            list: List containing phase numbers.
        """
        return self.ctf.get_phases_numbers()
    
    def phases_df(self):
        """
        Retrieves DataFrame containing phase counts and percentages.

        Returns:
            pandas.DataFrame: DataFrame containing phase counts and percentages.
        """
        df_phases = self.ctf.get_data()[0]
        phases_list = list(self.phases_names()["phase"])
        phases_list.insert(0, 'NaN')
        phase_names = {i: item for i, item in enumerate(phases_list)}
        df_phases['Phase'] = df_phases['Phase'].map(phase_names)
        phase_counts = df_phases['Phase'].value_counts().sort_index()
        total_count = phase_counts.sum()
        phase_percentages = (phase_counts / total_count) * 100
        phase_counts_sorted = phase_counts.sort_values(ascending=False)
        phase_percentages_sorted = phase_percentages.sort_values(ascending=False)
        result_df = pd.concat([phase_counts_sorted, phase_percentages_sorted], axis=1)
        result_df.columns = ['Count', 'Percentage']
        return result_df
    
    def get_ebsd_data_header(self):
        """
        Retrieves EBSD data and header information.

        Returns:
            tuple: Tuple containing EBSD data and header information.
        """
        data, header_data = self.ctf.get_data()
        return data, header_data
    
    def get_ebsd_data(self):
        """
        Retrieves EBSD data.

        Returns:
            pandas.DataFrame: DataFrame containing EBSD data.
        """
        data, header_data = self.ctf.get_data()
        return data
    
    def downsampleEBSD(self, factor=10, df=None):
        """
        Downsamples EBSD data.

        Parameters:
            factor (int, optional): Downsampling factor. Defaults to 10.
            df (pandas.DataFrame, optional): DataFrame containing EBSD data. Defaults to None.

        Returns:
            pandas.DataFrame: Downsampled DataFrame.
        """
        if df is None:
            df = self.get_ebsd_data()
        df = df[df.index % factor == 0]
        return df
    
    def phases(self, df_phases=None):
        """
        Retrieves phase percentages.

        Parameters:
            df_phases (pandas.DataFrame, optional): DataFrame containing phase information. Defaults to None.

        Returns:
            dict: Dictionary containing phase percentages.
        """
        if df_phases is None:
            df_phases = self.ctf.get_data()[0]
        phases_list = list(self.phases_names()["phase"])
        phases_list.insert(0, 'NaN')
        phase_counts = df_phases['Phase'].value_counts()
        phase_counts = phase_counts.sort_index()
        total_count = phase_counts.sum()
        phase_percentages = {phase: (count / total_count) * 100 for phase, count in zip(phases_list, phase_counts)}
        return phase_percentages
    
    def filterByPhase(self, phase_list, data=None):
        """
        Filters EBSD data by phase.

        Parameters:
            phase_list (list): List of phase names.
            data (pandas.DataFrame, optional): DataFrame containing EBSD data. Defaults to None.

        Returns:
            pandas.DataFrame: Filtered DataFrame.
        """
        if data is None:
            data = self.ctf.get_data()[0]
        indexes = []
        for phase in phase_list:
            phase_items = list(self.phases())
            index = phase_items.index(phase)
            indexes.append(index)
        for index in indexes:
            data = data[data['Phase'] != index]
        return data
    
    def getAnisotropyForEBSD(self, cij, euler_angles):
        tensor = Tensor()
        tensor_list = []
        for voigt in cij:
            tensor_list.append(tensor.voigt_to_tensor(voigt))

        rotated_tensor_list = []
        len_euler = []
        x = 0
        for euler_angle in euler_angles:
            len_euler.append(len(euler_angle))
            for i in range(len(euler_angle)):
                alpha = euler_angle.iloc[i]["Euler1"]
                beta = euler_angle.iloc[i]["Euler2"]
                gamma = euler_angle.iloc[i]["Euler3"]
                output = np.array(tensor.rotate_tensor(tensor_list[x], alpha, beta, gamma))
            #     print(output)
            #     print(f"Counter: {counter}")
            #     counter +=1
                rotated_tensor_list.append(output)
            x+= 1

        tensor_sum = np.sum(rotated_tensor_list, axis=0)
        tensor_sum= tensor_sum/sum(len_euler)

        return tensor_sum, len_euler
    

    def euler_to_quaternion(self, phi, phi1, phi2):
        """
        Convert Euler angles (Bunge angles) to quaternion.

        Parameters:
        phi : float
            First Euler angle (in radians).
        phi1 : float
            Second Euler angle (in radians).
        phi2 : float
            Third Euler angle (in radians).

        Returns:
        quaternion : numpy.ndarray
            Quaternion as [w, x, y, z].
        """
        r = Rotation.from_euler('ZXZ', [phi, phi1, phi2], degrees=True)
        return r.as_quat()



if __name__ == "__main__":
    # Example usage
    ctfobj = EBSD("ebsd.ctf")
    print(ctfobj.get_euler_angles(3))

    phi = 45
    phi1 = 30
    phi2 = 60 

    quaternion = euler_to_quaternion(phi, phi1, phi2)
    print("Quaternion:", quaternion)
