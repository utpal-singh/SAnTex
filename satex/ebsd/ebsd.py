# sage/ebsd/ebsd.py

from .ctf_parser import Ctf
import matplotlib.pyplot as plt 
from tabulate import tabulate

import numpy as np 
import pandas as pd 
from sklearn.neighbors import KNeighborsClassifier 
from joblib import Parallel, delayed 
from ..tensor import Tensor
from scipy.spatial.transform import Rotation

from .calcGrainBoundary import assign_to_grains_parallel
import plotly.express as px

from .melt import calcMelttensor
from .rotateEBSD import apply_custom_rotation_to_dataframe, apply_custom_rotation_to_dataframe_noxy

# from .ebsdrotation import apply_custom_rotation_to_dataframe as rotebsd
from .ebsdrotation import plot as plotrotebsd

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
        self.ctf = Ctf(self._filename) 
    
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
        

    def plot(self, data=None, rotation_angle=0, inside_plane=True, mirror=False):
        """
        Plots the EBSD map with colors based on phase. Allows for rotation and optional mirroring of the data.

        Parameters:
            data (pandas.DataFrame, optional): DataFrame containing EBSD data. If not provided, uses stored data.
            rotation_angle (int, optional): Angle by which to rotate the EBSD data (in degrees). Accepts 0, 90, 180, 270.
            inside_plane (bool, optional): If True, rotates the EBSD data inside the plane. If False, rotates outside the plane. Default is True.
            mirror (bool, optional): If True, mirrors the EBSD data horizontally before rotating. Default is False.

        Returns:
            None
        """
        if data is None:
            data, _ = self.ctf.get_data()

        # Mirror EBSD data if required
        if mirror:
            data['X'] = -data['X']  # Mirror horizontally

        # Rotate EBSD data
        if rotation_angle != 0:
            theta = np.radians(rotation_angle)
            if inside_plane:
                rotation_matrix = np.array([[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]])
            else:
                rotation_matrix = np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
            xy = np.column_stack((data['X'], data['Y']))
            xy_rotated = np.dot(xy, rotation_matrix)
            data['X'], data['Y'] = xy_rotated[:, 0], xy_rotated[:, 1]
        
        # Filter data based on phase
        df_new = data[data['Phase'].isin([0, 1, 2, 3])]
        
        # Define custom colormap for phases
        phases_colors = {0: 'black', 1: 'green', 2: 'yellow', 3:'blue'}
        
        # Plot
        fig, ax = plt.subplots(figsize=(16, 9))  # Adjust the figsize as per your requirement
        ax.scatter(df_new['X'], df_new['Y'], c=df_new['Phase'].map(phases_colors), s=1)  # Set s=1 for one point per pixel
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title('EBSD Map with Color Based on Phase')
        
        ax.set_aspect('equal')  # Set aspect ratio to 'equal' to ensure equal scaling along both axes

        plt.show()




    def get_index_of_phases(self, phases_list):
        phases = self.phases()
        phases_data = self.phases_data

        try:
            phases_list = ['Forsterite', 'Diopside']  
            phases_list = [phase.lower() for phase in phases_list]

            phase_indexes = [phase[0] for phase in phases_data if phase[1].lower() in phases_list]

            return phase_indexes

        except Exception as e:
            print("An error occurred:", e)


        user_phases_input = ['Forsterite', 'Diopside']

        phase_indexes = [phase[0] for phase in phases if phase[1] in user_phases_input]

        print(phase_indexes)


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
            str: Tabulated representation of phase percentages.
        """
        if df_phases is None:
            df_phases = self.ctf.get_data()[0]
        phases_list = list(self.phases_names()["phase"])
        phases_list.insert(0, 'NaN')
        phase_counts = df_phases['Phase'].value_counts()
        phase_counts = phase_counts.sort_index()
        total_count = phase_counts.sum()
        phase_percentages = {phase: (count / total_count) * 100 for phase, count in zip(phases_list, phase_counts)}
        
        data = [(index, phase, percentage) for index, (phase, percentage) in enumerate(phase_percentages.items())]
        
        headers = ["Index", "Phase", "Percentage"]

        self.phases_data = data
        
        table = tabulate(data, headers=headers, tablefmt="grid")
        return table


    
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
    
    def getAnisotropyForEBSD(self, cij, euler_angles, density, melt=0):
        if melt:
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
                    rotated_tensor_list.append(output)
                x+= 1

            density_averaged = (np.sum(np.multiply(density,len_euler)))/(np.sum(len_euler))

            tensor_sum = np.sum(rotated_tensor_list, axis=0)
            tensor_sum= tensor_sum/sum(len_euler)

            tensor_sum = (1 - melt*0.01)*tensor_sum + melt*0.01*calcMelttensor

            return tensor_sum, len_euler, density_averaged

        tensor = Tensor()
        tensor_list = []
        for voigt in cij:
            tensor_list.append(tensor.voigt_to_tensor(voigt))
            
        tensor = Tensor()
        rotated_tensor = []

        total_len_euler = []

        for j in (euler_angles):
            total_len_euler.append(len(j))

        a = np.array(total_len_euler)
        b = np.array(density)
        average_density = sum(b*a)/sum(a)

        for i in range(len(euler_angles)):
            for index, row in euler_angles[i].iterrows():
                euler1_value = row['Euler1']
                euler2_value = row['Euler2']
                euler3_value = row['Euler3']
            
                rotated_tensor.append(np.array((tensor.rotate_tensor(tensor.voigt_to_tensor(cij[i]), euler1_value, euler2_value, euler3_value))))

        average_tensor = np.mean(rotated_tensor, axis=0)

        return tensor.tensor_to_voigt(average_tensor), average_density
    

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
    
    def calcGrains(self, df, threshold, phase_names, downsampling_factor=20):
        """
        this is the df from df = ebsdfile.get_ebsd_data(), 
        threshold will be the threshold angle in degrees, phase_name will be a 
        list of phases present and that can be obtained from self.phases_names
        by default threshold is set to 10 degrees, but can be modified as threshold = 20 and
        so on and so forth
        """
        if df is None:
            df = self.get_ebsd_data()

        if phase_names is None:
            phase_names = self.phases_names()['phase'].tolist()
            phase_names.insert(0, "na")

        if threshold is None:
            threshold = 10

        # Filter out phase 0
        df = df[df['Phase'] != 0]

        if downsampling_factor:
            df = df.iloc[::downsampling_factor]
        else:
            df

        grain_indices = assign_to_grains_parallel(df[['Euler1', 'Euler2', 'Euler3']], threshold)
        df['Grain'] = pd.Series(grain_indices)
        for grain_idx, group in df.groupby('Grain'):
            phase_counts = group['Phase'].value_counts()
            dominant_phase = phase_counts.idxmax()
            print(f"Grain {grain_idx}: Dominant Phase - {phase_names[dominant_phase]}, Size - {len(group)}")

        self.plotGrains(df)
        return df
    
    def filterByGrainSize(self, df, phases_names, min_grain_size = 100):
        """
        df is the input from calcGrainsdf, with grain column appended
        """

        if df is None:
            df = self.calcGrains()

        if phases_names == None:
            phases_names = self.phases_names()['phase'].tolist()
            phases_names.insert(0, "na")        

        # Filter out grains below the minimum size threshold
        valid_grains = df['Grain'].value_counts()[df['Grain'].value_counts() >= min_grain_size].index

        # Filter the dataframe to include only rows corresponding to valid grains
        df_filtered = df[df['Grain'].isin(valid_grains)]

        # Display grains with their respective phase names
        for grain_idx, group in df_filtered.groupby('Grain'):
            phase_counts = group['Phase'].value_counts()
            dominant_phase = phase_counts.idxmax()
            print(f"Grain {grain_idx}: Dominant Phase - {phases_names[dominant_phase]}, Size - {len(group)}")

        self.plotGrains(df_filtered)

        return df_filtered
    
    def plotGrains(self, df):
        fig = px.scatter(df, x='X', y='Y', color='Phase', color_continuous_scale='viridis')

        # Customize marker size (2 pixel here)
        fig.update_traces(marker=dict(size=2))

        # Customize the layout
        fig.update_layout(
            title="2D Scatter Plot with Phase",
            xaxis_title="X",
            yaxis_title="Y",
            coloraxis_colorbar_title="Phase",
            showlegend=True
        )

        fig.show()

    def rotateEBSD(self, ebsd_df, angles, keepXY = True):

        if keepXY == False:
            return apply_custom_rotation_to_dataframe_noxy(ebsd_df, angles)
        else:
            return apply_custom_rotation_to_dataframe(ebsd_df, angles)
    

    def plot_rotate_ebsd(self, sample_ref = ["x2east", "zOutOfPlane"], ebsd_df = None, keepXY = False):
        if ebsd_df is None:
            ebsd_df = self.get_ebsd_data()

        if sample_ref == ["x2east", "zOutOfPlane"]:
            angle = (0, 0, 0)

        if sample_ref == ["x2west", "zOutOfPlane"]:
            angle = (0, 0, 180)

        if sample_ref == ["x2north", "zOutOfPlane"]:
            angle = (0, 0, 90)

        if sample_ref == ["x2south", "zOutOfPlane"]:
            angle = (0, 0, 270)

        if sample_ref == ["x2east", "zIntoPlane"]:
            angle = (0, 180, 0)

        if sample_ref == ["x2west", "zIntoPlane"]:
            angle = (0, 180, 180)

        if sample_ref == ["x2north", "zIntoPlane"]:
            angle = (0, 180, 90)

        if sample_ref == ["x2south", "zIntoPlane"]:
            angle = (0, 180, 270)

        if keepXY == False:
            rotated_ebsd_df = apply_custom_rotation_to_dataframe_noxy(ebsd_df, angle)
        else:
            rotated_ebsd_df = apply_custom_rotation_to_dataframe(ebsd_df, angle)

        return rotated_ebsd_df

        

        





if __name__ == "__main__":
    ctfobj = EBSD("ebsd.ctf")
    print(ctfobj.get_euler_angles(3))

    phi = 45  # First Euler angle
    phi1 = 30 # Second Euler angle
    phi2 = 60 # Third Euler angle

    quaternion = euler_to_quaternion(phi, phi1, phi2)
    print("Quaternion:", quaternion)