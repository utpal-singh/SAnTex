import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde

from orix.crystal_map import Phase
from orix.quaternion import Orientation, symmetry

from orix.vector import Miller, Vector3d
import plotly.graph_objects as go

import random

def pf(df, phase=1, crystal_symmetry='D2', random_val=True, uvw=[1, 0, 0], hemisphere = 'both', 
        axes_labels=["Xs", "Ys"], alpha = 0.01, figure = None, vector_labels = None, reproject=False, show_hemisphere_label = None,
        grid = None, grid_resolution = None, return_figure = None):
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
        - uvw: list of int, optional
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
    
    euler_df = df[df['Phase'] == phase][['Euler1', 'Euler2', 'Euler3']]

    euler_values = [tuple(x) for x in euler_df.to_numpy()]
    if random_val == True:
        random.shuffle(euler_values)
        euler_values = euler_values[:1250]
    else:
        euler_values = euler_values

    if crystal_symmetry in ["triclinic", "C1", "1", "-1", "1"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C1)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C1))
    elif crystal_symmetry in ["triclinic", "Ci", "-1", "-1", "1"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Ci)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.Ci))
    elif crystal_symmetry in ["monoclinic", "C2", "211", "2/m11", "211"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C2))
    elif crystal_symmetry in ["monoclinic", "Cs", "m11", "2/m11", "211"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Cs)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.Cs))
    elif crystal_symmetry in ["monoclinic", "C2h", "12/m1", "12/m1", "121"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C2h))
    elif crystal_symmetry in ["monoclinic", "C2", "121", "12/m1", "121"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C2))
    elif crystal_symmetry in ["monoclinic", "Cs", "1m1", "12/m1", "121"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Cs)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.Cs))
    elif crystal_symmetry in ["monoclinic", "C2h", "12/m1", "12/m1", "121"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C2h))
    elif crystal_symmetry in ["monoclinic", "C2", "112", "112/m", "112"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C2))
    elif crystal_symmetry in ["monoclinic", "Cs", "11m", "112/m", "112"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Cs)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.Cs))
    elif crystal_symmetry in ["monoclinic", "C2h", "112/m", "112/m", "112"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C2h))
    elif crystal_symmetry in ["orthorhombic", "D2", "222", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D2)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D2))
    elif crystal_symmetry in ["orthorhombic", "C2v", "2mm", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2v)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C2v))
    elif crystal_symmetry in ["orthorhombic", "C2v", "m2m", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2v)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C2v))
    elif crystal_symmetry in ["orthorhombic", "C2v", "mm2", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2v)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C2v))
    elif crystal_symmetry in ["orthorhombic", "D2h", "mmm", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D2h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D2h))
    elif crystal_symmetry in ["trigonal", "C3", "3", "-3", "3"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C3))
    elif crystal_symmetry in ["trigonal", "C3i", "-3", "-3", "3"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3i)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C3i))
    elif crystal_symmetry in ["trigonal", "D3", "321", "-3m1", "321"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D3))
    elif crystal_symmetry in ["trigonal", "C3v", "3m1", "-3m1", "321"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3v)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C3v))
    elif crystal_symmetry in ["trigonal", "D3d", "-3m1", "-3m1", "321"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3d)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D3d))
    elif crystal_symmetry in ["trigonal", "D3", "312", "-31m", "312"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D3))
    elif crystal_symmetry in ["trigonal", "C3v", "31m", "-31m", "312"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3v)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C3v))
    elif crystal_symmetry in ["trigonal", "D3d", "-31m", "-31m", "312"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3d)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D3d))
    elif crystal_symmetry in ["tetragonal", "C4", "4", "4/m", "4"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C4)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C4))
    elif crystal_symmetry in ["tetragonal", "S4", "-4", "4/m", "4"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.S4)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.S4))
    elif crystal_symmetry in ["tetragonal", "C4h", "4/m", "4/m", "4"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C4h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C4h))
    elif crystal_symmetry in ["tetragonal", "D4", "422", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D4)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D4))
    elif crystal_symmetry in ["tetragonal", "C4v", "4mm", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C4v)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C4v))
    elif crystal_symmetry in ["tetragonal", "D2d", "-42m", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D2d)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D2d))
    elif crystal_symmetry in ["tetragonal", "D2d", "-4m2", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D2d)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D2d))
    elif crystal_symmetry in ["tetragonal", "D4h", "4/mmm", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D4h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D4h))
    elif crystal_symmetry in ["hexagonal", "C6", "6", "6/m", "6"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C6)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C6))
    elif crystal_symmetry in ["hexagonal", "C3h", "-6", "6/m", "6"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C3h))
    elif crystal_symmetry in ["hexagonal", "C6h", "6/m", "6/m", "6"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C6h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C6h))
    elif crystal_symmetry in ["hexagonal", "D6", "622", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D6)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D6))
    elif crystal_symmetry in ["hexagonal", "C6v", "6mm", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C6v)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C6v))
    elif crystal_symmetry in ["hexagonal", "D3h", "-62m", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D3h))
    elif crystal_symmetry in ["hexagonal", "D3h", "-6m2", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D3h))
    elif crystal_symmetry in ["hexagonal", "D6h", "6/mmm", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D6h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D6h))
    elif crystal_symmetry in ["cubic", "T", "23", "m-3", "23"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.T)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.T))
    elif crystal_symmetry in ["cubic", "Th", "m-3", "m-3", "23"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Th)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.Th))
    elif crystal_symmetry in ["cubic", "O", "432", "m-3m", "432"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.O)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.O))
    elif crystal_symmetry in ["cubic", "Td", "-43m", "m-3m", "432"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Td)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.Td))
    elif crystal_symmetry in ["cubic", "Oh", "m-3m", "m-3m", "432"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Oh)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.Oh))

    uvw_ = uvw_.symmetrise(unique=False)

    uvw_ = uvw_.reshape(1, uvw_.size)
    orientations = orientations.reshape(*orientations.shape, 1)
    r_ = ~orientations * uvw_

    r_.scatter(hemisphere=hemisphere, axes_labels=axes_labels, alpha = alpha, figure = figure, vector_labels = vector_labels, reproject=reproject, show_hemisphere_label = show_hemisphere_label,
        grid = grid, grid_resolution = grid_resolution, return_figure = return_figure)
    return r_

def ipf(df, phase=1, vector_sample=[0, 0, 1], random_val=True,
            vector_title='Z', crystal_symmetry='D2'):
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
    
    from orix.quaternion import Orientation, symmetry
    from orix.vector import Miller, Vector3d
    import matplotlib.pyplot as plt
    from orix.crystal_map import Phase

    euler_df = df[df['Phase'] == phase][['Euler1', 'Euler2', 'Euler3']]

    euler_values = [tuple(x) for x in euler_df.to_numpy()]
    if random_val == True:
        random.shuffle(euler_values)
        euler_values = euler_values[:1250]
    else:
        euler_values = euler_values

    if crystal_symmetry in ["triclinic", "C1", "1", "-1", "1"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C1)
        symmetry_ = symmetry.C1
    elif crystal_symmetry in ["triclinic", "Ci", "-1", "-1", "1"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Ci)
        symmetry_ = symmetry.Ci
    elif crystal_symmetry in ["monoclinic", "C2", "211", "2/m11", "211"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2)
        symmetry_ = symmetry.C2
    elif crystal_symmetry in ["monoclinic", "Cs", "m11", "2/m11", "211"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Cs)
        symmetry_ = symmetry.Cs
    elif crystal_symmetry in ["monoclinic", "C2h", "12/m1", "12/m1", "121"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2h)
        symmetry_ =symmetry.C2h
    elif crystal_symmetry in ["monoclinic", "C2", "121", "12/m1", "121"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2)
        symmetry_ =symmetry.C2
    elif crystal_symmetry in ["monoclinic", "Cs", "1m1", "12/m1", "121"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Cs)
        symmetry_ =symmetry.Cs
    elif crystal_symmetry in ["monoclinic", "C2h", "12/m1", "12/m1", "121"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2h)
        symmetry_ =symmetry.C2h
    elif crystal_symmetry in ["monoclinic", "C2", "112", "112/m", "112"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2)
        symmetry_ =symmetry.C2
    elif crystal_symmetry in ["monoclinic", "Cs", "11m", "112/m", "112"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Cs)
        symmetry_ =symmetry.Cs
    elif crystal_symmetry in ["monoclinic", "C2h", "112/m", "112/m", "112"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2h)
        symmetry_ =symmetry.C2h
    elif crystal_symmetry in ["orthorhombic", "D2", "222", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D2)
        symmetry_ =symmetry.D2
    elif crystal_symmetry in ["orthorhombic", "C2v", "2mm", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2v)
        symmetry_ =symmetry.C2v
    elif crystal_symmetry in ["orthorhombic", "C2v", "m2m", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2v)
        symmetry_ =symmetry.C2v
    elif crystal_symmetry in ["orthorhombic", "C2v", "mm2", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2v)
        symmetry_ =symmetry.C2v
    elif crystal_symmetry in ["orthorhombic", "D2h", "mmm", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D2h)
        symmetry_ =symmetry.D2h
    elif crystal_symmetry in ["trigonal", "C3", "3", "-3", "3"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3)
        symmetry_ =symmetry.C3
    elif crystal_symmetry in ["trigonal", "C3i", "-3", "-3", "3"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3i)
        symmetry_ =symmetry.C3i
    elif crystal_symmetry in ["trigonal", "D3", "321", "-3m1", "321"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3)
        symmetry_ =symmetry.D3
    elif crystal_symmetry in ["trigonal", "C3v", "3m1", "-3m1", "321"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3v)
        symmetry_ =symmetry.C3v
    elif crystal_symmetry in ["trigonal", "D3d", "-3m1", "-3m1", "321"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3d)
        symmetry_ =symmetry.D3d
    elif crystal_symmetry in ["trigonal", "D3", "312", "-31m", "312"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3)
        symmetry_ =symmetry.D3
    elif crystal_symmetry in ["trigonal", "C3v", "31m", "-31m", "312"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3v)
        symmetry_ =symmetry.C3v
    elif crystal_symmetry in ["trigonal", "D3d", "-31m", "-31m", "312"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3d)
        symmetry_ =symmetry.D3d
    elif crystal_symmetry in ["tetragonal", "C4", "4", "4/m", "4"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C4)
        symmetry_ =symmetry.C4
    elif crystal_symmetry in ["tetragonal", "S4", "-4", "4/m", "4"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.S4)
        symmetry=symmetry.S4
    elif crystal_symmetry in ["tetragonal", "C4h", "4/m", "4/m", "4"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C4h)
        symmetry_ =symmetry.C4h
    elif crystal_symmetry in ["tetragonal", "D4", "422", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D4)
        symmetry_ =symmetry.D4
    elif crystal_symmetry in ["tetragonal", "C4v", "4mm", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C4v)
        symmetry_ =symmetry.C4v
    elif crystal_symmetry in ["tetragonal", "D2d", "-42m", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D2d)
        symmetry_ =symmetry.D2d
    elif crystal_symmetry in ["tetragonal", "D2d", "-4m2", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D2d)
        symmetry_ =symmetry.D2d
    elif crystal_symmetry in ["tetragonal", "D4h", "4/mmm", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D4h)
        symmetry_ =symmetry.D4h
    elif crystal_symmetry in ["hexagonal", "C6", "6", "6/m", "6"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C6)
        symmetry_ =symmetry.C6
    elif crystal_symmetry in ["hexagonal", "C3h", "-6", "6/m", "6"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3h)
        symmetry_ =symmetry.C3h
    elif crystal_symmetry in ["hexagonal", "C6h", "6/m", "6/m", "6"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C6h)
        symmetry_ =symmetry.C6h
    elif crystal_symmetry in ["hexagonal", "D6", "622", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D6)
        symmetry_ =symmetry.D6
    elif crystal_symmetry in ["hexagonal", "C6v", "6mm", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C6v)
        symmetry_ =symmetry.C6v
    elif crystal_symmetry in ["hexagonal", "D3h", "-62m", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3h)
        symmetry_ =symmetry.D3h
    elif crystal_symmetry in ["hexagonal", "D3h", "-6m2", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3h)
        symmetry_ =symmetry.D3h
    elif crystal_symmetry in ["hexagonal", "D6h", "6/mmm", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D6h)
        symmetry_ =symmetry.D6h
    elif crystal_symmetry in ["cubic", "T", "23", "m-3", "23"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.T)
        symmetry_ =symmetry.T
    elif crystal_symmetry in ["cubic", "Th", "m-3", "m-3", "23"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Th)
        symmetry_ =symmetry.Th
    elif crystal_symmetry in ["cubic", "O", "432", "m-3m", "432"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.O)
        symmetry_ =symmetry.O
    elif crystal_symmetry in ["cubic", "Td", "-43m", "m-3m", "432"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Td)
        symmetry_ =symmetry.Td
    elif crystal_symmetry in ["cubic", "Oh", "m-3m", "m-3m", "432"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Oh)
        symmetry_ =symmetry.Oh

    orientations = orientations.reshape(*orientations.shape, 1)
    # print(orientations)
    vec_sample = Vector3d(vector_sample)
    # print(vec_sample)
    vec_crystal = orientations * vec_sample
    subplot_kw = dict(projection='ipf', symmetry=symmetry_)
    fig = plt.figure(figsize=(9, 8))

    # ax0 = fig.add_subplot(221, direction=vec_sample, **subplot_kw)
    # ax0.scatter(orientations, alpha=0.05)
    # _ = ax0.set_title(f"Phase{phase}, {vector_title}")

    ax = fig.add_subplot(223, direction=vec_sample, **subplot_kw)
    ax.pole_density_function(vec_crystal)
    _ = ax.set_title(f"Phase{phase}, {vector_title}")


def pdf(df, phase=1, crystal_symmetry='D2', random_val=True, uvw=[0, 1, 0], hemisphere = 'both', sigma = 4,
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
        - uvw: list of int, optional
            The Miller indices representing the crystallographic plane to calculate ODF for (default is [1, 0, 0]).
        - axes_labels: list of str or None, optional
            Labels for the X and Y axes of the ODF plot (default is None).
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
    
    euler_df = df[df['Phase'] == phase][['Euler1', 'Euler2', 'Euler3']]

    euler_values = [tuple(x) for x in euler_df.to_numpy()]
    if random_val == True:
        random.shuffle(euler_values)
        euler_values = euler_values[:1250]
    else:
        euler_values = euler_values

    if crystal_symmetry in ["triclinic", "C1", "1", "-1", "1"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C1)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C1))
    elif crystal_symmetry in ["triclinic", "Ci", "-1", "-1", "1"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Ci)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.Ci))
    elif crystal_symmetry in ["monoclinic", "C2", "211", "2/m11", "211"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C2))
    elif crystal_symmetry in ["monoclinic", "Cs", "m11", "2/m11", "211"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Cs)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.Cs))
    elif crystal_symmetry in ["monoclinic", "C2h", "12/m1", "12/m1", "121"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C2h))
    elif crystal_symmetry in ["monoclinic", "C2", "121", "12/m1", "121"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C2))
    elif crystal_symmetry in ["monoclinic", "Cs", "1m1", "12/m1", "121"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Cs)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.Cs))
    elif crystal_symmetry in ["monoclinic", "C2h", "12/m1", "12/m1", "121"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C2h))
    elif crystal_symmetry in ["monoclinic", "C2", "112", "112/m", "112"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C2))
    elif crystal_symmetry in ["monoclinic", "Cs", "11m", "112/m", "112"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Cs)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.Cs))
    elif crystal_symmetry in ["monoclinic", "C2h", "112/m", "112/m", "112"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C2h))
    elif crystal_symmetry in ["orthorhombic", "D2", "222", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D2)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D2))
    elif crystal_symmetry in ["orthorhombic", "C2v", "2mm", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2v)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C2v))
    elif crystal_symmetry in ["orthorhombic", "C2v", "m2m", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2v)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C2v))
    elif crystal_symmetry in ["orthorhombic", "C2v", "mm2", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2v)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C2v))
    elif crystal_symmetry in ["orthorhombic", "D2h", "mmm", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D2h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D2h))
    elif crystal_symmetry in ["trigonal", "C3", "3", "-3", "3"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C3))
    elif crystal_symmetry in ["trigonal", "C3i", "-3", "-3", "3"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3i)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C3i))
    elif crystal_symmetry in ["trigonal", "D3", "321", "-3m1", "321"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D3))
    elif crystal_symmetry in ["trigonal", "C3v", "3m1", "-3m1", "321"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3v)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C3v))
    elif crystal_symmetry in ["trigonal", "D3d", "-3m1", "-3m1", "321"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3d)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D3d))
    elif crystal_symmetry in ["trigonal", "D3", "312", "-31m", "312"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D3))
    elif crystal_symmetry in ["trigonal", "C3v", "31m", "-31m", "312"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3v)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C3v))
    elif crystal_symmetry in ["trigonal", "D3d", "-31m", "-31m", "312"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3d)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D3d))
    elif crystal_symmetry in ["tetragonal", "C4", "4", "4/m", "4"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C4)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C4))
    elif crystal_symmetry in ["tetragonal", "S4", "-4", "4/m", "4"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.S4)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.S4))
    elif crystal_symmetry in ["tetragonal", "C4h", "4/m", "4/m", "4"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C4h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C4h))
    elif crystal_symmetry in ["tetragonal", "D4", "422", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D4)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D4))
    elif crystal_symmetry in ["tetragonal", "C4v", "4mm", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C4v)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C4v))
    elif crystal_symmetry in ["tetragonal", "D2d", "-42m", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D2d)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D2d))
    elif crystal_symmetry in ["tetragonal", "D2d", "-4m2", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D2d)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D2d))
    elif crystal_symmetry in ["tetragonal", "D4h", "4/mmm", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D4h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D4h))
    elif crystal_symmetry in ["hexagonal", "C6", "6", "6/m", "6"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C6)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C6))
    elif crystal_symmetry in ["hexagonal", "C3h", "-6", "6/m", "6"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C3h))
    elif crystal_symmetry in ["hexagonal", "C6h", "6/m", "6/m", "6"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C6h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C6h))
    elif crystal_symmetry in ["hexagonal", "D6", "622", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D6)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D6))
    elif crystal_symmetry in ["hexagonal", "C6v", "6mm", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C6v)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.C6v))
    elif crystal_symmetry in ["hexagonal", "D3h", "-62m", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D3h))
    elif crystal_symmetry in ["hexagonal", "D3h", "-6m2", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D3h))
    elif crystal_symmetry in ["hexagonal", "D6h", "6/mmm", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D6h)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.D6h))
    elif crystal_symmetry in ["cubic", "T", "23", "m-3", "23"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.T)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.T))
    elif crystal_symmetry in ["cubic", "Th", "m-3", "m-3", "23"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Th)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.Th))
    elif crystal_symmetry in ["cubic", "O", "432", "m-3m", "432"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.O)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.O))
    elif crystal_symmetry in ["cubic", "Td", "-43m", "m-3m", "432"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Td)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.Td))
    elif crystal_symmetry in ["cubic", "Oh", "m-3m", "m-3m", "432"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Oh)
        uvw_ = Miller(uvw = uvw, phase=Phase(point_group=symmetry.Oh))

    uvw_ = uvw_.symmetrise(unique=False)

    uvw_ = uvw_.reshape(1, uvw_.size)
    orientations = orientations.reshape(*orientations.shape, 1)
    r_ = ~orientations * uvw_
    r_.pole_density_function(sigma=sigma, hemisphere=hemisphere, axes_labels=axes_labels, figure = figure, show_hemisphere_label = show_hemisphere_label,
        grid = grid, grid_resolution = grid_resolution, return_figure = return_figure, log = log, colorbar=colorbar, weights = weights)


def calculate_odf(euler_angles):
    """
    Calculate and plot an ODF from Bunge Euler angles.
    
    Parameters:
    - euler_angles: A Nx3 numpy array of Euler angles (phi1, Phi, phi2) in degrees.
    
    """
    euler_angles_rad = np.radians(euler_angles)
    
    # using a Gaussian KDE to estimate the distribution
    kde = gaussian_kde(euler_angles_rad.T)
    
    phi1, Phi, phi2 = np.mgrid[0:np.pi:100j, 0:np.pi/2:100j, 0:2*np.pi:100j]
    positions = np.vstack([phi1.ravel(), Phi.ravel(), phi2.ravel()])
    
    density = np.reshape(kde(positions).T, phi1.shape)
    
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(phi1, Phi, phi2, c=density, cmap='viridis')
    ax.set_xlabel('Phi1')
    ax.set_ylabel('Phi')
    ax.set_zlabel('Phi2')
    ax.set_title('ODF Representation')
    plt.show()


def plot_orientations_3d(df):

    euler_angles = df[['Euler1', 'Euler2', 'Euler3']].to_numpy()


    phi1, Phi, phi2 = euler_angles[:, 0], euler_angles[:, 1], euler_angles[:, 2]
    
    fig = go.Figure(data=[go.Scatter3d(
        x=phi1,
        y=Phi,
        z=phi2,
        mode='markers',
        marker=dict(
            size=5,
            color='rgba(255, 0, 0, 0.8)',  
        )
    )])
    
    fig.update_layout(
        title='3D Plot of Bunge-Euler Orientations',
        scene=dict(
            xaxis_title='Phi1 (Degrees)',
            yaxis_title='Phi (Degrees)',
            zaxis_title='Phi2 (Degrees)'
        ),
        margin=dict(r=0, b=0, l=0, t=30) 
    )
    
    fig.show()


def ipf_colorkey(df, phase=1, crystal_symmetry='D2'):
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
    
    from orix.quaternion import Orientation, symmetry
    from orix.vector import Miller, Vector3d
    import matplotlib.pyplot as plt
    from orix.crystal_map import Phase
    from orix import plot, sampling

    euler_df = df[df['Phase'] == phase][['Euler1', 'Euler2', 'Euler3']]

    euler_values = [tuple(x) for x in euler_df.to_numpy()]

    if crystal_symmetry in ["triclinic", "C1", "1", "-1", "1"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C1)
        point_group=symmetry.C1
    elif crystal_symmetry in ["triclinic", "Ci", "-1", "-1", "1"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Ci)
        point_group=symmetry.Ci
    elif crystal_symmetry in ["monoclinic", "C2", "211", "2/m11", "211"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2)
    elif crystal_symmetry in ["monoclinic", "Cs", "m11", "2/m11", "211"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Cs)
        point_group=symmetry.Cs
    elif crystal_symmetry in ["monoclinic", "C2h", "12/m1", "12/m1", "121"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2h)
        point_group=symmetry.C2h
    elif crystal_symmetry in ["monoclinic", "C2", "121", "12/m1", "121"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2)
        point_group=symmetry.C2
    elif crystal_symmetry in ["monoclinic", "Cs", "1m1", "12/m1", "121"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Cs)
        point_group=symmetry.Cs
    elif crystal_symmetry in ["monoclinic", "C2h", "12/m1", "12/m1", "121"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2h)
        point_group=symmetry.C2h
    elif crystal_symmetry in ["monoclinic", "C2", "112", "112/m", "112"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2)
        point_group=symmetry.C2
    elif crystal_symmetry in ["monoclinic", "Cs", "11m", "112/m", "112"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Cs)
        point_group=symmetry.Cs
    elif crystal_symmetry in ["monoclinic", "C2h", "112/m", "112/m", "112"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2h)
        point_group=symmetry.C2h
    elif crystal_symmetry in ["orthorhombic", "D2", "222", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D2)
        point_group=symmetry.D2
    elif crystal_symmetry in ["orthorhombic", "C2v", "2mm", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2v)
        point_group=symmetry.C2v
    elif crystal_symmetry in ["orthorhombic", "C2v", "m2m", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2v)
        point_group=symmetry.C2v
    elif crystal_symmetry in ["orthorhombic", "C2v", "mm2", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C2v)
        point_group=symmetry.C2v
    elif crystal_symmetry in ["orthorhombic", "D2h", "mmm", "mmm", "222"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D2h)
        point_group=symmetry.D2h
    elif crystal_symmetry in ["trigonal", "C3", "3", "-3", "3"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3)
        point_group=symmetry.C3
    elif crystal_symmetry in ["trigonal", "C3i", "-3", "-3", "3"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3i)
        point_group=symmetry.C3i
    elif crystal_symmetry in ["trigonal", "D3", "321", "-3m1", "321"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3)
        point_group=symmetry.D3
    elif crystal_symmetry in ["trigonal", "C3v", "3m1", "-3m1", "321"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3v)
        point_group=symmetry.C3v
    elif crystal_symmetry in ["trigonal", "D3d", "-3m1", "-3m1", "321"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3d)
        point_group=symmetry.D3d
    elif crystal_symmetry in ["trigonal", "D3", "312", "-31m", "312"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3)
        point_group=symmetry.D3
    elif crystal_symmetry in ["trigonal", "C3v", "31m", "-31m", "312"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3v)
        point_group=symmetry.C3v
    elif crystal_symmetry in ["trigonal", "D3d", "-31m", "-31m", "312"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3d)
        point_group=symmetry.D3d
    elif crystal_symmetry in ["tetragonal", "C4", "4", "4/m", "4"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C4)
        point_group=symmetry.C4
    elif crystal_symmetry in ["tetragonal", "S4", "-4", "4/m", "4"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.S4)
        point_group=symmetry.S4
    elif crystal_symmetry in ["tetragonal", "C4h", "4/m", "4/m", "4"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C4h)
        point_group=symmetry.C4h
    elif crystal_symmetry in ["tetragonal", "D4", "422", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D4)
        point_group=symmetry.D4
    elif crystal_symmetry in ["tetragonal", "C4v", "4mm", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C4v)
        point_group=symmetry.C4v
    elif crystal_symmetry in ["tetragonal", "D2d", "-42m", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D2d)
        point_group=symmetry.D2d
    elif crystal_symmetry in ["tetragonal", "D2d", "-4m2", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D2d)
        point_group=symmetry.D2d
    elif crystal_symmetry in ["tetragonal", "D4h", "4/mmm", "4/mmm", "422"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D4h)
        point_group=symmetry.D4h
    elif crystal_symmetry in ["hexagonal", "C6", "6", "6/m", "6"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C6)
        point_group=symmetry.C6
    elif crystal_symmetry in ["hexagonal", "C3h", "-6", "6/m", "6"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C3h)
        point_group=symmetry.C3h
    elif crystal_symmetry in ["hexagonal", "C6h", "6/m", "6/m", "6"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C6h)
        point_group=symmetry.C6h
    elif crystal_symmetry in ["hexagonal", "D6", "622", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D6)
        point_group=symmetry.D6
    elif crystal_symmetry in ["hexagonal", "C6v", "6mm", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.C6v)
        point_group=symmetry.C6v
    elif crystal_symmetry in ["hexagonal", "D3h", "-62m", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3h)
        point_group=symmetry.D3h
    elif crystal_symmetry in ["hexagonal", "D3h", "-6m2", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D3h)
        point_group=symmetry.D3h
    elif crystal_symmetry in ["hexagonal", "D6h", "6/mmm", "6/mmm", "622"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.D6h)
        point_group=symmetry.D6h
    elif crystal_symmetry in ["cubic", "T", "23", "m-3", "23"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.T)
        point_group=symmetry.T
    elif crystal_symmetry in ["cubic", "Th", "m-3", "m-3", "23"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Th)
        point_group=symmetry.Th
    elif crystal_symmetry in ["cubic", "O", "432", "m-3m", "432"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.O)
        point_group=symmetry.O
    elif crystal_symmetry in ["cubic", "Td", "-43m", "m-3m", "432"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Td)
        point_group=symmetry.Td
    elif crystal_symmetry in ["cubic", "Oh", "m-3m", "m-3m", "432"]:
        orientations = Orientation.from_euler(np.deg2rad(euler_values), symmetry.Oh)
        point_group=symmetry.Oh

    ipfkey = plot.IPFColorKeyTSL(point_group)
    orientations.symmetry = ipfkey.symmetry
    rgb_z = ipfkey.orientation2color(orientations)
    orientations.scatter("ipf", c=rgb_z, direction=ipfkey.direction)

