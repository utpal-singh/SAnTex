# sage/ebsd/ebsd.py

from ._ctf_parser import Ctf
from matplotlib import pyplot as plt
from tabulate import tabulate
import os

import numpy as np 
import pandas as pd 
from ..tensor import Tensor
from scipy.spatial.transform import Rotation

from .calcGrainBoundary import assign_to_grains_parallel
import plotly.express as px

from .melt_tensor import calcMelttensor
from .rotateEBSD import apply_custom_rotation_to_dataframe, apply_custom_rotation_to_dataframe_noxy

# from .ebsdrotation import apply_custom_rotation_to_dataframe as rotebsd
from .odf import ipf, odf, pdf

class EBSD:
    """
    Class for handling Electron Backscatter Diffraction (EBSD) data.
    """

    def __init__(self, filename) -> None:
        """ 
        Initializes the EBSD object.

        Parameters:
            filename (str): The filename of the EBSD data file.

        Raises:
            FileNotFoundError: If the specified file does not exist.
        """
        self._filename = filename
        if not os.path.isfile(self._filename):
                    raise FileNotFoundError(f"The file '{self._filename}' was not found.")
        
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
        

    def plot(self, data=None, rotation_angle=0, inside_plane=True, mirror=False, save_image=False, image_filename=None, 
            dpi=300, cmap='viridis', legend_location="upper right"):
        """
        Plots the EBSD map with colors based on phase. Allows for rotation and optional mirroring of the data.

        Parameters:
            data (pandas.DataFrame, optional): DataFrame containing EBSD data. If not provided, uses stored data.
            rotation_angle (int, optional): Angle by which to rotate the EBSD data (in degrees). Accepts 0, 90, 180, 270.
            inside_plane (bool, optional): If True, rotates the EBSD data inside the plane. If False, rotates outside the plane. Default is True.
            mirror (bool, optional): If True, mirrors the EBSD data horizontally before rotating. Default is False.
            save_image (bool, optional): If True, saves the plot as an image. Default is False.
            image_filename (str, optional): Filename to use when saving the image. Required if save_image is True.
            dpi (int, optional): Dots per inch for the saved image. Default is 300.
            cmap (str, optional): Colormap to use for plotting. Default is 'viridis'.
            legend_location (str, optional): Location of the legend. Options are 'upper right', 'upper left', 'lower right', 'lower left'. Default is 'upper right'.

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
        df_new = data

        # Plot
        fig, ax = plt.subplots(figsize=(16, 9))  # Adjust the figsize as per your requirement
        scatter = ax.scatter(df_new['X'], df_new['Y'], c=df_new['Phase'], cmap=cmap, s=1)
        ax.set_xlabel('X')
        ax.set_ylabel('Y')
        ax.set_title('EBSD Map with Color Based on Phase')

        ax.set_aspect('equal')  # Set aspect ratio to 'equal' to ensure equal scaling along both axes

        # Create legend for phase colors
        handles, _ = scatter.legend_elements()
        legend = ax.legend(handles, [f'Phase {phase}' for phase in sorted(df_new['Phase'].unique())], title="Phases", loc=legend_location)
        ax.add_artist(legend)

        if save_image and image_filename:
            plt.savefig(image_filename, dpi=dpi)  # Save the plot as an image with user-specified filename and DPI

        plt.show()


    def get_index_of_phases(self, phases_list):
        """
        Get the index of phases based on their names.

        Parameters:
            phases_list (list): A list of phase names.

        Returns:
            list: A list containing the indexes of the specified phases.
        """
        phases = self.phases()
        phases_data = self.phases_data

        try:
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

        Parameters:
            phase (str): Phase name.

        Returns:
            pandas.dataframe: euler angles of phase by name
        """
        phase_index = self.get_index_of_phases(phase)
        return self.get_euler_angles(phase_index)

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
        
        return data

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
            index = next(index for index, mineral, value in phase_items if mineral == phase)
            indexes.append(index)
        for index in indexes:
            data = data[data['Phase'] != index]
        return data
    
    def getAnisotropyForEBSD(self, cij, euler_angles, density, melt=0):
        """
        Calculates average density and tensor for EBSD data

        Parameters:
        cij: array, stiffness tensor of phase
        euler_angles: list, euler angles of phases in list
        density: float, density of phases
        melt: float, fraction of melt available

        Returns:
        - average_tensor: array, average tensor calculated
        - density: float, average density
        """
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

            return tensor_sum, density_averaged

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
    
    def calcGrains(self, df, phase_names, threshold = 10, downsampling_factor=20):
        """
        Calculated the grains in the ebsd data

        Parameters:
        - df: EBSD dataframe, this is the df we get from from df = ebsdfile.get_ebsd_data(), 
        - threshold: float, threshold will be the threshold angle in degrees, threshold mean orientation angle, by default threshold is set to 10 degrees, but can be modified as threshold = 20 and so on and so forth
        - phase_names: list, phase_name will be a list of phases present and that can be obtained from self.phases_names
        - downsampling factor: int, specify downsampling factor for faster computation, default is 20

        Returns:
        - df: EBSD dataframe with grains calculated
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
        Filters ebsd by grain size

        Parameters:
        - df: EBSD dataframe, which is the input from calcGrainsdf, with grain column appended
        - min_grain_size: int, desired minimum grain size, default is 100

        Returns:
        - df_filtered: Filtered EBSD dataframe
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
    
    def plotGrains(self, df, color_continuous_scale='viridis', save_name=None, dpi=None):
        """
        Create a 2D scatter plot with phase coloring.

        Parameters:
        - df (DataFrame): The EBSD DataFrame containing the data to plot.
        - color_continuous_scale (str): The color scale to use for continuous colors, default is viridis.
        - save_name (str): The name to use when saving the figure as a PNG file. If None, the figure will be displayed but not saved.
        - dpi (int): The DPI (dots per inch) to use when saving the figure. If None, the default DPI will be used.

        Returns:
        - None
        """
        fig = px.scatter(df, x='X', y='Y', color='Phase', color_continuous_scale=color_continuous_scale)
        fig.update_traces(marker=dict(size=2))
        fig.update_layout(
            title="2D Scatter Plot with Phase",
            xaxis_title="X",
            yaxis_title="Y",
            coloraxis_colorbar_title="Phase",
            showlegend=True
        )
        if save_name:
            if not save_name.endswith('.png'):
                save_name += '.png'  # Ensure the file has a .png extension
            if dpi:
                fig.write_image(save_name, format='png', scale=dpi)
            else:
                fig.write_image(save_name, format='png')
            print(f"Figure saved as {os.path.abspath(save_name)}")
        else:
            fig.show()

    def rotateEBSD(self, ebsd_df, angles, keepXY = True):

        """
        Rotate EBSD with a certain angles

        Parameters:
        - ebsd_df: ebsd pandas dataframe
        - angles: angles to be rotated (phi1, phi2, phi3), eg: [20, 30, 40]
        - keepXY: bool, Specify if you want to keep X and Y constant. Default is True

        Returns:
        - Rotated EBSD dataframe 
        
        """

        if keepXY == False:
            return apply_custom_rotation_to_dataframe_noxy(ebsd_df, angles)
        else:
            return apply_custom_rotation_to_dataframe(ebsd_df, angles)
    

    def plot_rotate_ebsd(self, sample_ref = ["x2east", "zOutOfPlane"], ebsd_df = None, keepXY = False):
        """
        Plots the rotated ebsd and returns the rotated ebsd dataframe

        Parameters:
        - sample_ref: sample reference, eg: ["x2east", "zOutOfPlane"], ["x2west", "zOutOfPlane"], ["x2north", "zOutOfPlane"], ["x2south", "zOutOfPlane"], ["x2east", "zIntoPlane"], default is ["x2east", "zOutOfPlane"]
        - ebsd_df: the ebsd pandas dataframe
        - keepXY: bool, Specify if you want to keep X and Y constant. Default is False

        Returns:
        - rotated_ebsd_df : rotated ebsd dataframe
        """
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
    
    def filterMAD(self, df, value=0.7):
        """
        Filters the ebsd dataframe with MAD greater than a certain threshold

        Parameters:
        - df: pandas EBSD dataframe
        - value: value of the threshold for MAD (default is 0.7)

        Returns:
        - df: EBSD filtered dataframe
        """
        return df[df['MAD']<value]
    
    def filterByPhaseNumberList(self, phase_list, df=None):
        """
        Filters the EBSD dataframe given a list of phases

        Parameters:
        - df: pandas EBSD dataframe
        - phase_list: list of the index of phases, eg: phase_list = [4, 5, 6, 7]
        """
        if df is None:
            df = self.ctf.get_data()[0]

        df_filtered = df[~df['Phase'].isin(phase_list)]
        return df_filtered
    
    def odf(self, df, phase=1, crystal_symmetry='D2', random_val=True, miller=[1, 0, 0], hemisphere = 'both', 
            axes_labels=["Xs", "Ys"], alpha = 0.01, figure = None, vector_labels = None, reproject=False, show_hemisphere_label = None,
            grid = None, grid_resolution = None, return_figure = None):
        """
            Calculate the Orientation Distribution Function (ODF) for a given EBSD dataset and optionally visualize it.

            Parameters:
            - df: pandas EBSD DataFrame
                The DataFrame containing EBSD (Electron Backscatter Diffraction) data.
            - phase: int, optional
                The phase number for which ODF needs to be calculated (default is 1).
            - crystal_symmetry: str, optional
                The crystal symmetry of the material, e.g., 'D2' for cubic symmetry (default is 'D2').
            - random_val: bool, optional
                If True, randomly generate orientation values for incomplete data (default is True).
            - miller: list of int, optional
                The Miller indices representing the crystallographic plane to calculate ODF for (default is [1, 0, 0]).
            - projection: str, optional
                The type of projection to use for visualization (default is 'stereographic').
            - figure: matplotlib Figure or None, optional
                If provided, the ODF plot will be added to this Figure; otherwise, a new Figure will be created.
            - axes_labels: list of str or None, optional
                Labels for the X and Y axes of the ODF plot (default is None).
            - vector_labels: list of str or None, optional
                Labels for the vectors in the ODF plot (default is None, which results in auto-generated labels).
            - hemisphere: str, upper or lower or both, optional
                Specify the hemisphere for ODF calculation and visualization (default is None).
            - reproject: bool, optional
                If True, reproject data into a specified coordinate system (default is False).
            - show_hemisphere_label: bool or None, optional
                Specify whether to show the hemisphere label on the plot (default is None).
            - grid: bool or None, optional
                Specify whether to display a grid on the ODF plot (default is None).
            - grid_resolution: tuple of float or None, optional
                Resolution of the grid (default is None).
            - figure_kwargs: dict or None, optional
                Additional keyword arguments for configuring the matplotlib Figure (default is None).
            - reproject_scatter_kwargs: dict or None, optional
                Additional keyword arguments for configuring scatter plot during reprojection (default is None).
            - text_kwargs: dict or None, optional
                Additional keyword arguments for configuring text elements in the plot (default is None).
            - return_figure: bool, optional
                If True, return the matplotlib Figure object along with the ODF data (default is False).

            Returns:
            - odf_result: pandas DataFrame
                DataFrame containing the calculated Orientation Distribution Function (ODF) values.
            - figure: matplotlib Figure or None
                If return_figure is True, the Figure object containing the ODF plot; otherwise, None.
            """
        odf(df=df, phase=phase, crystal_symmetry=crystal_symmetry, random_val=random_val, miller=miller, hemisphere = hemisphere, 
            axes_labels=axes_labels, alpha = alpha, figure = figure, vector_labels = vector_labels, reproject=reproject, show_hemisphere_label = show_hemisphere_label,
            grid = grid, grid_resolution = grid_resolution, return_figure = return_figure)

    def pdf(self, df, phase=1, crystal_symmetry='D2', random_val=True, miller=[0, 1, 0], hemisphere = 'both', sigma = 4,
        axes_labels=["Xs", "Ys"], figure = None, show_hemisphere_label = None,
        grid = None, grid_resolution = None, return_figure = None, log = False, colorbar=True, weights = None):

        """
            Calculate the Orientation Distribution Function (ODF) for a given EBSD dataset and optionally visualize it.

            Parameters:
            - df: pandas DataFrame
                The DataFrame containing EBSD (Electron Backscatter Diffraction) data.
            - phase: int, optional
                The phase number for which ODF needs to be calculated (default is 1).
            - crystal_symmetry: str, optional
                The crystal symmetry of the material, e.g., 'D2' for cubic symmetry (default is 'D2').
            - random_val: bool, optional
                If True, randomly generate orientation values for incomplete data (default is True).
            - miller: list of int, optional
                The Miller indices representing the crystallographic plane to calculate ODF for (default is [1, 0, 0]).
            - projection: str, optional
                The type of projection to use for visualization (default is 'stereographic').
            - figure: matplotlib Figure or None, optional
                If provided, the ODF plot will be added to this Figure; otherwise, a new Figure will be created.
            - axes_labels: list of str or None, optional
                Labels for the X and Y axes of the ODF plot (default is None).
            - vector_labels: list of str or None, optional
                Labels for the vectors in the ODF plot (default is None, which results in auto-generated labels).
            - hemisphere: str, upper or lower or both, optional
                Specify the hemisphere for ODF calculation and visualization (default is None).
            - reproject: bool, optional
                If True, reproject data into a specified coordinate system (default is False).
            - show_hemisphere_label: bool or None, optional
                Specify whether to show the hemisphere label on the plot (default is None).
            - grid: bool or None, optional
                Specify whether to display a grid on the ODF plot (default is None).
            - grid_resolution: tuple of float or None, optional
                Resolution of the grid (default is None).
            - figure_kwargs: dict or None, optional
                Additional keyword arguments for configuring the matplotlib Figure (default is None).
            - reproject_scatter_kwargs: dict or None, optional
                Additional keyword arguments for configuring scatter plot during reprojection (default is None).
            - text_kwargs: dict or None, optional
                Additional keyword arguments for configuring text elements in the plot (default is None).
            - return_figure: bool, optional
                If True, return the matplotlib Figure object along with the ODF data (default is False).

            Returns:
            - odf_result: pandas DataFrame
                DataFrame containing the calculated Orientation Distribution Function (ODF) values.
            - figure: matplotlib Figure or None
                If return_figure is True, the Figure object containing the ODF plot; otherwise, None.
            """

        pdf(df = df, phase=phase, crystal_symmetry=crystal_symmetry, random_val=random_val, miller=miller, sigma=sigma, hemisphere=hemisphere, axes_labels=axes_labels, figure = figure, show_hemisphere_label = show_hemisphere_label,
        grid = grid, grid_resolution = grid_resolution, return_figure = return_figure, log = log, colorbar=colorbar, weights = weights)

    def ipf(self, df, phase=1, vector_sample=[0, 0, 1], random_val=True,
            vector_title='Z', projection='ipf', crystal_symmetry='D2'):
        
        """
            Generate and visualize an Inverse Pole Figure (IPF) from EBSD data.

            Parameters:
            - df: pandas DataFrame
                The DataFrame containing EBSD (Electron Backscatter Diffraction) data.
            - phase: int, optional
                The phase number for which IPF needs to be generated (default is 1).
            - vector_sample: list of int, optional
                The sample vector used for IPF generation (default is [0, 0, 1]).
            - random_val: bool, optional
                If True, randomly select orientation values for IPF generation (default is True).
            - vector_title: str, optional
                The title for the vector used in IPF visualization (default is 'Z').
            - projection: str, optional
                The type of projection for IPF visualization (default is 'ipf').
            - crystal_symmetry: str, optional
                The crystal symmetry of the material (default is 'D2').

            Returns:
            - fig: matplotlib Figure
                The generated matplotlib Figure object containing the IPF visualization.
        """

        self.pdf(df = df, phase=phase, random_val=random_val, crystal_symmetry=crystal_symmetry)
        ipf(df = df, phase=phase, vector_sample=vector_sample, random_val=random_val,
            vector_title=vector_title, projection=projection, crystal_symmetry=crystal_symmetry)
