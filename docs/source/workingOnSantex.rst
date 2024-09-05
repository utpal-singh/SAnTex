Working with SAnTex
========================

In a Python working environment, import SAnTex's modules and Python libraries:

.. code-block:: python

    from santex.ebsd import EBSD
    from santex.tensor import Tensor
    from santex.material import Material
    from santex.isotropy import Isotropy
    import pandas as pd


Load the EBSD file:

.. code-block:: python

    ebsdfile = EBSD("ebsd.ctf")


This loads the `ebsd.ctf` file in the Python object `ebsdfile` of class `EBSD`, from which users are able to call methods for further processing, cleaning of the ebsd data, and ultimately calculate elastic properties.


This loads the ebsd.ctf file in the python object ebsdfile of class EBSD, from which
users are able to call methods for further processing, cleaning of the ebsd data, and
ultimately calculate elastic properties.

*Class:* In Python, a class is a blueprint for creating objects. It defines the structure
and behavior of objects of that type. Think of it as a template or a prototype that
encapsulates data (attributes) and methods (functions) that operate on that data.

*Object:* An object, also known as an instance, is a unique entity that is created based
on the structure defined by its class. When an object is instantiated, it inherits the
attributes and methods defined in its class, but can also have its own unique state
and behavior.

**Note:** The `ebsd.ctf` should be located in the directory this notebook is being run from. But you may provide a relative path like ``data/ebsd.ctf`` or an absolute import such as ``/Users/myname/Documents/ebsd/ebsd.ctf``.


Available methods in the EBSD class
==========

To list the phases present in the ebsd file:


.. code-block:: python

    phases_available = ebsdfile.phases()
    print(phases_available)

Load the ebsd file into a pandas dataframe ``df``, via invoking the ``get_ebsd_data()``
method from the EBSD class: *Note:* this dataframe ``df`` is critical as it holds
information of all the euler angles and the ebsd data.


.. code-block:: python

    df = ebsdfile.get_ebsd_data()

To look into the contents of the `df`, enter the following directive:

.. code-block:: python

    print(df)

To look into the data header of the EBSD file, following can be used:

.. code-block:: python

    ebsd_data_header = ebsdfile.get_ebsd_data_header()
    print(ebsd_data_header)

To look at the index of phase, phases, and relative abundence, follwoing can be entered:

.. code-block:: python

    phases = ebsdfile.phases()
    print(phases)

To get the euler angles of a certain phase, lets say phase = 1 in this case is Forsterite, enter the following:

.. code-block:: python

    forsteriteEulerAngles = ebsdfile.get_euler_angles(phase = 1, data = df)
    print(forsteriteEulerAngles)

To plot this ebsd file, following can be enterred:

.. code-block:: python

    ebsdfile.plot(df)

To save this image, following can be used:

.. code-block:: python

    ebsdfile.plot(df, cmap = "viridis", save_image= True, image_filename= "ebsd.jpg", legend_location="lower right")


Plotting Conventions
===================

Following are the keywords to orient sample reference plane and to store in a new
dataframe, and plot it

.. code-block:: bash

    1. sample_ref = ["x2east", "zOutOfPlane"]
    2. sample_ref = ["x2west", "zOutOfPlane"]
    3. sample_ref = ["x2north", "zOutOfPlane"]
    4. sample_ref = ["x2south", "zOutOfPlane"]
    5. sample_ref = ["x2east", "zIntoPlane"]
    6. sample_ref = ["x2west", "zIntoPlane"]
    7. sample_ref = ["x2north", "zIntoPlane"]
    8. sample_ref = ["x2south", "zIntoPlane"]

To rotate the ebsd dataframe to x2east and zIntoPlane, do the following:

.. code-block:: python

    rotated_df = ebsdfile.plot_rotate_ebsd(sample_ref = ["x2east",
    "zIntoPlane"], ebsd_df = df)

This rotates the ebsdfile in the specified convention, and saves the rotated dataframe in user-defined variable ``rotated_df``

To plot this rotated_df, user can enter the following command:

.. code-block:: python

    rotated_df = ebsdfile.plot_rotate_ebsd(sample_ref = ["x2east",
    "zIntoPlane"], ebsd_df = df)


**Note:** ebsdfile object created initaially from ``ebsdfile = EBSD('ebsd.ctf')`` can be reused for methods of ``EBSD`` object. For example, to plot the rotated dataframe, following can be entered:

.. code-block:: python

    ebsdfile.plot(rotated_df)


To rotate the EBSD data to match the SEM orientation in any custom angles, a user can enter the angles in Bunze ZXZ format as:

.. code-block:: python

    angles = (180, 0, 0)
    updatedebsd = ebsdfile.rotateEBSD(rotated_df, angles)

To view the current orientations, a user can always look at the ebsd plot using:

.. code-block:: python

    ebsdfile.plot(updatedebsd)

Get the phases names using the following:

.. code-block:: python

    phases_names = ebsdfile.phases_names()
    print(phases_names)

Cleaning EBSD data
===========

In this example, to clean the EBSD dataset we successively remove grains with large mean angular
deviation (MAD) of 0.8, reconstruct grains with grain boundaries misorientation ≥ 10
degrees, and remove grains smaller than 7 pixels.

Remove the phases which are mis-indexed (enter the endex of the phases):
--------------

.. code-block:: python

    df = ebsdfile.filterByPhaseNumberList(df = updatedebsd, phase_list = [4, 5, 6, 7])

Remove pixel data with MAD higher than a user specified value
---------

To remove pixels with mean angular deviation e.g. (MAD) > 0.8, and store the
cleaned dataset into a new dataframe called filtered_df:

.. code-block:: python

    filtered_df = ebsdfile.filterMAD(df, 0.8)


Reconstruct grains with boundaries whose misorientation exceeds a minimum of 10 degrees:
---------------------

To reconstruct grains with misorientation ≥ 10 degrees and store the cleaner dataset
into a new dataframe called df_grain_boundary:

.. code-block:: python

    phases_names = ebsdfile.phases_names()
    phases_names = phases_names['phase'].tolist()
    phases_names.insert(0, "na")

    df_grain_boundary = ebsdfile.calcGrains(df = filtered_df, threshold = 10, phase_names=phases_names, downsampling_factor=10)

**Note:** Here the downsampling_factor is applied just for speedy handson!

Remove small grains
---------------

User can remove grains smaller than e.g. 7 pixels, and store the cleaner dataset into a
new dataframe called filtered_df_grain_boundary:

.. code-block:: python

    filtered_df_grain_boundary = ebsdfile.filterByGrainSize(df_grain_boundary, phases_names, min_grain_size=7)

Compare original and clean datasets
----------------


.. code-block:: python

    ebsdfile.plot()
    ebsdfile.plot(df = filtered_df_grain_boundary)


Calculating Anisotropy from EBSD file:
==============

To prepare the dataframe, the user can instanciate Material class as follows:

.. code-block:: python

    material_instance = Material()

The densities are within the santex package registry, and can be accessed via the following directive. The pressure in GPa and Temperature can be entered as:

.. code-block:: python

    rho_Fo = material_instance.load_density("Forsterite", 2, 1500)
    rho_diop = material_instance.load_density("Diopside", 2, 1500)
    rho_ens = material_instance.load_density("Enstatite", 2, 1500)

The registry for elastic tensors within santex can be accessed via the following commands. The pressure and temperacture
conditions can be entered via the keyword PRESSURE and TEMP

.. code-block:: python

    cij_Fo = material_instance.voigthighPT('Forsterite', PRESSURE = 2, TEMP = 1500)
    cij_ens = material_instance.voigthighPT('Enstatite', PRESSURE = 2, TEMP = 1500)
    cij_diop = material_instance.voigthighPT('Diopside', PRESSURE = 2, TEMP = 1500)

For preparing dataframes for seismic anisotropy, the elasticstiffness tensors and densities needs to be passed as a list, as:

.. code-block:: python

    cij = [cij_Fo, cij_ens, cij_diop]
    density = [rho_Fo, rho_ens, rho_diop]

The euler angles for the phases should also be bundled in list as:

.. code-block:: python

    forsterite = ebsdfile.get_euler_angles(phase = 1, data=filtered_df_grain_boundary)
    enstatite = ebsdfile.get_euler_angles(phase = 2, data=filtered_df_grain_boundary)
    diopside = ebsdfile.get_euler_angles(phase = 3, data=filtered_df_grain_boundary)
    euler_angles = [forsterite, enstatite, diopside]

The anisotropy can then be calculated as following:

.. code-block:: python

    average_tensor, average_density = ebsdfile.getAnisotropyForEBSD(cij, euler_angles, density)
    anis = Anisotropy(average_tensor*10**9, average_density)

To look at the plots for the seismic velocities, following command can be entered:

.. code-block:: python

    anis.plot()

To look at the anisotropy values, maxvp, minvp, maxvs1, minvs1, etc.., folllowing command can be used:

.. code-block:: python

    anis.anisotropy_values()

This gives values such as:

.. code-block:: bash

    Max Vp:  8449.545069720665
    Min Vp:  7834.1684457727015
    Max Vs1:  4843.772824149944
    Min Vs1:  4577.275239228743
    Max Vs2:  4661.180158160457
    Min Vs2:  4555.315881762446
    Max vs anisotropy percent:  6.110766862825683
    Min vs anisotropy percent:  0.020842429078130564
    P wave anisotropy percent:  7.558185340984474
    S1 Wave anisotropy percent:  5.657493372889703
    S2 Wave anisotropy percent:  2.297278183366879
    Velocity difference:  287.209223367252
    Vp/Vs1 ratio:  3.4012758430117644


Ordinary Distribution Function (ODF)
=========

The orientation distribution function (ODF) is a function on the orientation space that associates to each orientation g the volume percentage of crystals in a polycrystaline specimen that are in this specific orientation

ODF can be calculated and plotted as:

.. code-block:: python

    ebsdfile.odf(df = filtered_df_grain_boundary, phase=1, crystal_symmetry='D2',
            random_val=True,
            miller=[1, 0, 0],
            hemisphere='both',
            axes_labels=['Xs', 'Ys'],
            alpha=0.01,
            figure=None,
            vector_labels=None,
            reproject=False,
            show_hemisphere_label=None,
            grid=None,
            grid_resolution=None,
            return_figure=None)


Pole Figures
===========

Pole figures are two dimensional representations of orientations. To illustrate this we define a random orientation with trigonal crystal symmetry

The pole figures can be calculated and plotted as:

.. code-block:: python

    ebsdfile.pdf(df=filtered_df_grain_boundary,
        phase=1,
        crystal_symmetry='D2',
        random_val=True,
        miller=[0, 1, 0],
        hemisphere='both',
        sigma=4,
        axes_labels=['Xs', 'Ys'],
        figure=None,
        show_hemisphere_label=None,
        grid=None,
        grid_resolution=None,
        return_figure=None,
        log=False,
        colorbar=True,
        weights=None,)

Following is the conversion tables which can be found from https://mtex-toolbox.github.io/HomepageOld/files/doc/symmetry_index.html . The crystal symmetry above
should be defined as following convention. The user can enter the crystal symmetry as in either Schoenflies, or International, or Laue class, or Rotational axes convention.

.. code-block:: python

    id  crystal system  Schoen-  Inter-    Laue     Rotational
                        flies    national  class    axes
    1   triclinic       C1       1         -1       1
    2   triclinic       Ci       -1        -1       1
    3   monoclinic      C2       211       2/m11    211
    4   monoclinic      Cs       m11       2/m11    211
    5   monoclinic      C2h      2/m11     2/m11    211
    6   monoclinic      C2       121       12/m1    121
    7   monoclinic      Cs       1m1       12/m1    121
    8   monoclinic      C2h      12/m1     12/m1    121
    9   monoclinic      C2       112       112/m    112
    10  monoclinic      Cs       11m       112/m    112
    11  monoclinic      C2h      112/m     112/m    112
    12  orthorhombic    D2       222       mmm      222
    13  orthorhombic    C2v      2mm       mmm      222
    14  orthorhombic    C2v      m2m       mmm      222
    15  orthorhombic    C2v      mm2       mmm      222
    16  orthorhombic    D2h      mmm       mmm      222
    17  trigonal        C3       3         -3       3
    18  trigonal        C3i      -3        -3       3
    19  trigonal        D3       321       -3m1     321
    20  trigonal        C3v      3m1       -3m1     321
    21  trigonal        D3d      -3m1      -3m1     321
    22  trigonal        D3       312       -31m     312
    23  trigonal        C3v      31m       -31m     312
    24  trigonal        D3d      -31m      -31m     312
    25  tetragonal      C4       4         4/m      4
    26  tetragonal      S4       -4        4/m      4
    27  tetragonal      C4h      4/m       4/m      4
    28  tetragonal      D4       422       4/mmm    422
    29  tetragonal      C4v      4mm       4/mmm    422
    30  tetragonal      D2d      -42m      4/mmm    422
    31  tetragonal      D2d      -4m2      4/mmm    422
    32  tetragonal      D4h      4/mmm     4/mmm    422
    33  hexagonal       C6       6         6/m      6
    34  hexagonal       C3h      -6        6/m      6
    35  hexagonal       C6h      6/m       6/m      6
    36  hexagonal       D6       622       6/mmm    622
    37  hexagonal       C6v      6mm       6/mmm    622
    38  hexagonal       D3h      -62m      6/mmm    622
    39  hexagonal       D3h      -6m2      6/mmm    622
    40  hexagonal       D6h      6/mmm     6/mmm    622
    41  cubic           T        23        m-3      23
    42  cubic           Th       m-3       m-3      23
    43  cubic           O        432       m-3m     432
    44  cubic           Td       -43m      m-3m     432
    45  cubic           Oh       m-3m      m-3m     432
    46  icosahedral     I        532       -5-32/m  532
    47  icosahedral     Ih       -5-32/m   -5-32/m  532


Inverse Pole Figure
==========

For an orientation distribution function (ODF) :math:`f : SO(3) \rightarrow \mathbb{R}`, the inverse pole density function :math:`P_{\vec{r}}` with respect to a fixed specimen direction :math:`\vec{r}` is a spherical function defined as the integral:

.. math::

    P_{\vec{r}}(\vec{h}) = \int_{\vec{g}\cdot\vec{h}=\vec{r}} f(\vec{g}) \, d\vec{g}

The pole density function :math:`P_{\vec{r}}(\vec{h})` evaluated at a crystal direction :math:`\vec{h}` can be interpreted as the volume percentage of crystals with the crystal lattice planes :math:`\vec{h}` being normal to the specimen direction :math:`\vec{r}`.

In order to illustrate the concept of inverse pole figures, let's calculate the ipf and plot:


.. code-block:: python

    ebsdfile.ipf(df=filtered_df_grain_boundary,
        phase=1,
        vector_sample=[0, 0, 1],
        random_val=True,
        vector_title='Z',
        projection='ipf',
        crystal_symmetry='D2',)