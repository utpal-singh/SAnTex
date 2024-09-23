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


Material analysis
===================

The Material class from the Sage library is used for defining and working with materials calculations.

Material Module can be loaded as:

.. code-block:: python

    import numpy as np
    from tabulate import tabulate
    from santex.material import Material

The available material within the santex registry can be viewed as:

.. code-block:: python

    material_instance = Material()
    phases_info = material_instance.availablePhases()
    print("Available Phases:")
    print(phases_info)

This loads a list of materials present with their crystal systems and primary phase information as shown below:

.. code-block:: bash

    Available Phases:
    +---------------------------------------+---------------------+-----------------------------------+
    |                 Phase                 |   Crystal System    |           Primary Phase           |
    +---------------------------------------+---------------------+-----------------------------------+
    |           Almandine-pyrope            |        Cubic        |              Garnet               |
    |               Grossular               |        Cubic        |              Garnet               |
    |               Majorite                |        Cubic        |              Garnet               |
    |                Pyrope                 |        Cubic        |              Garnet               |
    |              a_quartz_1               | Hexagonal/ Trigonal |              Quartz               |
    |              a_quartz_2               | Hexagonal/ Trigonal |              Quartz               |
    |              a_Quartz_3               | Hexagonal/ Trigonal |              Quartz               |
    |              a_quartz_4               | Hexagonal/ Trigonal |              Quartz               |
    |             a_quartz_696C             | Hexagonal/ Trigonal |              Quartz               |
    |             a_quartz_700C             | Hexagonal/ Trigonal |              Quartz               |
    |               Calcite_1               | Hexagonal/ Trigonal |              Calcite              |
    |               Calcite_2               | Hexagonal/ Trigonal |              Calcite              |
    |              Forsterite               |    Orthorhombic     |              Olivine              |
    |               Fayalite                |    Orthorhombic     |              Olivine              |
    |               Lawsonite               |    Orthorhombic     |              Olivine              |
    |        Orthoenstatite (MgSiO3)        |    Orthorhombic     |              Olivine              |
    |        Orthoenstatite (MgSiO3)        |    Orthorhombic     |              Olivine              |
    |               Enstatite               |    Orthorhombic     |              Olivine              |
    |       Bronzite (Mg0.8Fe0.2SiO3)       |    Orthorhombic     |              Olivine              |
    |         Ferrosilite (FeSiO3)          |    Orthorhombic     |              Olivine              |
    |                Biotite                |     Monoclinic      | Phyllosilicates and clay minerals |
    |               Muscovite               |     Monoclinic      | Phyllosilicates and clay minerals |
    |              Phlogopite               |     Monoclinic      | Phyllosilicates and clay minerals |
    |            Illite-smectite            |     Monoclinic      | Phyllosilicates and clay minerals |
    |                Dickite                |     Monoclinic      | Phyllosilicates and clay minerals |
    |                Augite                 |     Monoclinic      |          Clinopyroxenes           |
    |               Diopside                |     Monoclinic      |          Clinopyroxenes           |
    |            Chrome-diopside            |     Monoclinic      |          Clinopyroxenes           |
    |                Jadeite                |     Monoclinic      |          Clinopyroxenes           |
    |               Omphacite               |     Monoclinic      |          Clinopyroxenes           |
    |                Coesite                |     Monoclinic      |          Clinopyroxenes           |
    |       Amphobole #1 Richterite1        |     Monoclinic      |             Amphibole             |
    |       Amphobole #2 Kataphorite1       |     Monoclinic      |             Amphibole             |
    |  Amphobole #3 Taramite-Tschermakite1  |     Monoclinic      |             Amphibole             |
    | Amphobole #4 Hornblende-Tschermakite1 |     Monoclinic      |             Amphibole             |
    |        Amphobole #5 Tremolite1        |     Monoclinic      |             Amphibole             |
    |         Amphobole #6 Edenite1         |     Monoclinic      |             Amphibole             |
    |         Amphobole #7 Edenite1         |     Monoclinic      |             Amphibole             |
    |        Amphobole #8 Pargasite1        |     Monoclinic      |             Amphibole             |
    |        Amphobole #9 Pargasite1        |     Monoclinic      |             Amphibole             |
    |            Hornblende (#1)            |     Monoclinic      |             Amphibole             |
    |            Hornblende (#2)            |     Monoclinic      |             Amphibole             |
    |              Glaucophane              |     Monoclinic      |             Amphibole             |
    |          Sanidine (Or83Ab15)          |     Monoclinic      |         Alkali feldspars          |
    |          Sanidine (Or89Ab11)          |     Monoclinic      |         Alkali feldspars          |
    |         Orthoclase (Or93Ab7)          |     Monoclinic      |         Alkali feldspars          |
    |           Albite (Or0Ab100)           |      Triclinic      |       Plagioclase feldspar        |
    |             An0 (Albite)              |      Triclinic      |       Plagioclase feldspar        |
    |           An25 (Oligoclase)           |      Triclinic      |       Plagioclase feldspar        |
    |            An37 (Andesine)            |      Triclinic      |       Plagioclase feldspar        |
    |            An48 (Andesine)            |      Triclinic      |       Plagioclase feldspar        |
    |          An60 (Labradorite)           |      Triclinic      |       Plagioclase feldspar        |
    |          An67 (Labradorite)           |      Triclinic      |       Plagioclase feldspar        |
    |           An78 (Bytownite)            |      Triclinic      |       Plagioclase feldspar        |
    |           An96 (Anorthite)            |      Triclinic      |       Plagioclase feldspar        |
    |               Kaolinite               |      Triclinic      |           Clay minerals           |
    |                Nacrite                |      Triclinic      |           Clay minerals           |
    +---------------------------------------+---------------------+-----------------------------------+


Hooke's Law
===========

Hooke's Law describes the behavior of certain materials when subjected to a stretching or compressing force. Hooke's law can be expressed in terms of the elastic stiffness tensor and the strain tensor, as:

.. math::
   :label: eq:hookeslaw

   \sigma_{ij} = C_{ijkl} \epsilon_{kl}

where :math:`\sigma_{ij}` and :math:`\epsilon_{kl}` are the components of the stress and strain tensors, respectively, while :math:`C_{ijkl}` are the components of the elastic stiffness tensor. In this form, Hooke's law is more general and can account for the anisotropy and directionality of the material’s elastic properties.

The pressure and temperature dependence of elastic constants is mainly linear but can include non-linear effects that can be approximated up to second-order terms using a Taylor series expansion, as outlined below:

.. math::
   :label: eq:elasticity

   C_{ijkl}(p, T) = C_{ijkl}(0, 0) + \left. \frac{\partial C_{ijkl}}{\partial p} \right|_{(0,0)} p + \left. \frac{\partial C{ijkl}}{\partial T} \right|_{(0,0)} T + \mathcal{O}(p^2, T^2)

In the current version of SAnTex, melt is considered as an isotropic phase with homogenous distribution within an anisotropic host rock, e.g., [@lee_modeling_2017:2017].

.. math::
   :label: eq:elasticity_fmelt

   C_{ijkl}(p, T) = (1-f)\left(C_{ijkl}(0, 0) + \left. \frac{\partial C_{ijkl}}{\partial p} \right|_{(0,0)} p + \left. \frac{\partial C{ijkl}}{\partial T} \right|_{(0,0)} T + \mathcal{O}(p^2, T^2)\right) + f\text{melt}(C_{\text{melt}}(p, T))

The fraction of melt, :math:`f`, can be controlled by the user. :math:`C_{\text{melt}}` is the stiffness tensor of the melt. The approach currently incorporated in SAnTex overlooks the complex behavior of melt, including its viscosity, flow dynamics, and interaction with neighboring minerals, which can influence the overall anisotropic properties of the system. Future developments of SAnTex will aim to include more functionalities towards the calculation of melt-induced anisotropy.

Modal Mineral Composition of a rock and anisotropy
----------------

within Santex, a user can parse a modal mineral composition for a rock and calculate its anisotropy as follows:

.. code-block:: python

    from santex.material import Material
    from santex.anisotropy import Anisotropy

    material = Material()
    rock = ["Forsterite", "Diopside", "Enstatite"]
    fraction = [0.6, 0.25, 0.15]
    average_tensor, average_density = material.modalRock(rock, fraction, 2, 1000)
    anisotropy = Anisotropy(average_tensor*10**9, average_density)
    values = anisotropy.anisotropy_values()

This returns a dictionary containing the anisotropy percent containing ``maxvp``, ``minvp``, ``minvs1``, ``minvs2``, ``maxvs1``, ``maxvs2``, ``max_vs_anisotropy_percent``, ``min_vs_anisotropy_percent``, ``p_wave_anisotropy_percent``,
``s1_wave_anisotropy_percent``, ``s2_wave_anisotropy_percent``, ``maxdvs``, and ``AVpVs1``

To plot the seismic velocities:

.. code-block::python

    anisotropy.plot()


Isotropic Velocities analysis
----------------

Isotropy module can be loaded as follows:

.. code-block:: python

    from santex.isotropy import Isotropy
    isotropy = Isotropy()

To check the available phases, a user can invoke ``get_available_phases()`` within the Isotropy class as

.. code-block:: python

    isotropy.get_available_phases()

This gives a list of phases which are available in santex registry as is shown in below example

.. code-block:: bash

    ######################Available Phases######################
    Material id: aqz       Material name: Alpha-Quartz
    ######################Available Phases######################
    Material id: bqz       Material name: Beta-Quartz
    ######################Available Phases######################
    Material id: coe       Material name: Coesite
    ######################Available Phases######################
    Material id: hAb       Material name: High-T Albite
    ######################Available Phases######################
    Material id: lAb       Material name: Low-T Albite
    ######################Available Phases######################
    Material id: an       Material name: Anorthite
    ######################Available Phases######################
    Material id: or       Material name: Orthoclase
    ######################Available Phases######################
    Material id: san       Material name: Sanidine
    ######################Available Phases######################
    Material id: alm       Material name: Almandine
    ######################Available Phases######################
    Material id: gr       Material name: Grossular
    ######################Available Phases######################
    Material id: py       Material name: Pyrope
    ######################Available Phases######################
    Material id: fo       Material name: Forsterite
    ######################Available Phases######################
    Material id: fa       Material name: Fayalite
    ######################Available Phases######################
    Material id: en       Material name: Enstatite
    ######################Available Phases######################
    Material id: fs       Material name: Ferrosillite
    ######################Available Phases######################
    Material id: mgts       Material name: Mg-Tschermak
    ######################Available Phases######################
    Material id: di       Material name: Diopside
    ######################Available Phases######################
    Material id: hed       Material name: Hedenbergite
    ######################Available Phases######################
    Material id: jd       Material name: Jadeite
    ######################Available Phases######################
    Material id: ac       Material name: Acmite
    ######################Available Phases######################
    Material id: cats       Material name: Ca-Tschermak
    ######################Available Phases######################
    Material id: gl       Material name: Glaucophane
    ######################Available Phases######################
    Material id: fgl       Material name: Ferroglaucophane
    ######################Available Phases######################
    Material id: tr       Material name: Tremolite
    ######################Available Phases######################
    Material id: fact       Material name: Ferroactinolite
    ######################Available Phases######################
    Material id: ts       Material name: Tschermakite
    ######################Available Phases######################
    Material id: parg       Material name: Pargasite
    ######################Available Phases######################
    Material id: hb       Material name: Hornblende
    ######################Available Phases######################
    Material id: anth       Material name: Anthophyllite
    ######################Available Phases######################
    Material id: phl       Material name: Phlogopite
    ######################Available Phases######################
    Material id: ann       Material name: Annite
    ######################Available Phases######################
    Material id: mu       Material name: Muscovite
    ######################Available Phases######################
    Material id: cel       Material name: Celadonite
    ######################Available Phases######################
    Material id: ta       Material name: Talc
    ######################Available Phases######################
    Material id: clin       Material name: Clinochlore
    ######################Available Phases######################
    Material id: daph       Material name: Daphnite
    ######################Available Phases######################
    Material id: atg       Material name: Antigorite
    ######################Available Phases######################
    Material id: zo       Material name: Zoisite
    ######################Available Phases######################
    Material id: cz       Material name: Clinozoisite
    ######################Available Phases######################
    Material id: ep       Material name: Epidote
    ######################Available Phases######################
    Material id: law       Material name: Lawsonite
    ######################Available Phases######################
    Material id: pre       Material name: Prehnite
    ######################Available Phases######################
    Material id: pump       Material name: Pumpellyte
    ######################Available Phases######################
    Material id: lmt       Material name: Laumontite
    ######################Available Phases######################
    Material id: wrk       Material name: Wairakite
    ######################Available Phases######################
    Material id: br       Material name: Brucite
    ######################Available Phases######################
    Material id: chum       Material name: Clinohumite
    ######################Available Phases######################
    Material id: phA       Material name: Phase A
    ######################Available Phases######################
    Material id: sill       Material name: Sillimanite
    ######################Available Phases######################
    Material id: ky       Material name: Kyanite
    ######################Available Phases######################
    Material id: sp       Material name: Mg-Spinel
    ######################Available Phases######################
    Material id: herc       Material name: hercynite
    ######################Available Phases######################
    Material id: mt       Material name: Magnetite
    ######################Available Phases######################
    Material id: ilm       Material name: Ilmenite
    ######################Available Phases######################
    Material id: rut       Material name: Rutile
    ######################Available Phases######################
    Material id: ttn       Material name: Titanite
    ######################Available Phases######################
    Material id: crd       Material name: Cordierite
    ######################Available Phases######################
    Material id: scap       Material name: Scapolite
    ######################Available Phases######################
    Material id: cc       Material name: Calcite
    ######################Available Phases######################
    Material id: arag       Material name: Aragonite
    ######################Available Phases######################
    Material id: mag       Material name: Magnesite

Following input parameters are inbuilt in ythe santex registry which initiates the isotropic properties calculation. 

.. code-block:: python

    rho0: initial density
    ao: coefficient of thermal expansion
    akt0: isothermal bulk modulus, which is a measure of a material's resistance to compression under uniform pressure
    dkdp: pressure derivative of the bulk modulus, indicating how the bulk modulus changes with pressure
    amu0: shear modulus of the mineral. The shear modulus measures a material's resistance to deformation by shear stress
    dmudp: pressure derivative of the shear modulus, indicating how the shear modulus changes with pressure
    gam: gamma, first thermodynamic Gruinesen parameter
    grun: second Gruneisen parameter, which is a measure of how a material's volume changes with temperature
    delt: Debye temperature, which is a measure of the average vibrational energy of atoms in a solid.

User can call ``get_phase_constants()`` method as follows:

.. code-block:: python

    isotropy.get_phase_constants("Forsterite")

This returns a dictionary which looks like

.. code-block:: python

        {'id': 'fo',
    'name': 'Forsterite',
    'rho0': 3222.0,
    'ao': 6.13e-05,
    'akt0': 127300000000.0,
    'dkdp': 4.2,
    'amu0': 81600000000.0,
    'dmudp': 1.6,
    'gam': 5.19,
    'grun': 1.29,
    'delt': 5.5}

We can get the following quantities at any given temperature and pressure for a material after calling the method ``calculate_seismic_properties()``

1. density: density of material at any given pressure and temperature
2. aks: bulk modulus, The bulk modulus indicates how much a material will compress under pressure.
3. amu: Shear Modulus, The shear modulus is essential for understanding a material's response to shear stress
4. vp: P-wave velocity at any given pressure and temperature
5. vs: swave velocity at any given pressure and temperature
6. vbulk: Bulk sound velocity, The velocity of sound waves traveling through a material
7. akt: Isothermal bulk modulus, Similar to the bulk modulus, but specifically refers to the resistance to compression under constant 8. temperature conditions

The statement that can be used as:

.. code-block:: python
    density, aks, amu, vp, vs, vbulk, akt = isotropy.calculate_seismic_properties('Forsterite', temperature=2000, pressure=2, return_vp_vs_vbulk=True, return_aktout=True)


Tensor analysis
-------------

Stiffness tensors of phases are a 3*3*3*3 size, however, due to its symmetricity, it can also be written as a voigt notation of 6*6 matrix. A tensor can be rotated
in Bunge Euler (ZXZ) convention

To initaially instanciate the Tensor Class, following needs to be inputted in python:

.. code-block:: python

    import numpy as np
    from santex.tensor import Tensor

    tensor = Tensor()

Let's load some atiffness matrix values for forsterite written in voigt notation. The values can be seen in voigt matrix format

.. code-block:: python

    cij_forsterite = np.array([[320.5,  68.15,  71.6,     0,     0,     0],
    [ 68.15,  196.5,  76.8,     0,     0,     0],
    [  71.6,   76.8, 233.5,     0,     0,     0],
    [   0,      0,      0,     64,     0,     0],
    [   0,      0,      0,      0,    77,     0],
    [   0,      0,      0,      0,     0,  78.7]])

To convert this voigt notation to tensor notation, a user can call ``voigt_to_tensor()`` method within Tensor class as follows:

.. code-block:: python
    cijkl_forsterite = tensor.voigt_to_tensor(cij_forsterite)

This gives a 3*3*3*3 array with the converted voigt in full tensor notation as seen

.. code-block:: python

    array([[[[320.5 ,   0.  ,   0.  ],
         [  0.  ,  68.15,   0.  ],
         [  0.  ,   0.  ,  71.6 ]],

        [[  0.  ,  78.7 ,   0.  ],
         [ 78.7 ,   0.  ,   0.  ],
         [  0.  ,   0.  ,   0.  ]],

        [[  0.  ,   0.  ,  77.  ],
         [  0.  ,   0.  ,   0.  ],
         [ 77.  ,   0.  ,   0.  ]]],


       [[[  0.  ,  78.7 ,   0.  ],
         [ 78.7 ,   0.  ,   0.  ],
         [  0.  ,   0.  ,   0.  ]],

        [[ 68.15,   0.  ,   0.  ],
         [  0.  , 196.5 ,   0.  ],
         [  0.  ,   0.  ,  76.8 ]],

        [[  0.  ,   0.  ,   0.  ],
         [  0.  ,   0.  ,  64.  ],
         [  0.  ,  64.  ,   0.  ]]],


       [[[  0.  ,   0.  ,  77.  ],
         [  0.  ,   0.  ,   0.  ],
         [ 77.  ,   0.  ,   0.  ]],

        [[  0.  ,   0.  ,   0.  ],
         [  0.  ,   0.  ,  64.  ],
         [  0.  ,  64.  ,   0.  ]],

        [[ 71.6 ,   0.  ,   0.  ],
         [  0.  ,  76.8 ,   0.  ],
         [  0.  ,   0.  , 233.5 ]]]])

The converted tensor can be converted back to voigt notation as:

.. code-block:: python
    cij_forsterite = tensor.tensor_to_voigt(cijkl_forsterite)

Which then returns the 6*6 voigt notation of the tensor

.. code-block:: python

    array([[320.5 ,  68.15,  71.6 ,   0.  ,   0.  ,   0.  ],
       [ 68.15, 196.5 ,  76.8 ,   0.  ,   0.  ,   0.  ],
       [ 71.6 ,  76.8 , 233.5 ,   0.  ,   0.  ,   0.  ],
       [  0.  ,   0.  ,   0.  ,  64.  ,   0.  ,   0.  ],
       [  0.  ,   0.  ,   0.  ,   0.  ,  77.  ,   0.  ],
       [  0.  ,   0.  ,   0.  ,   0.  ,   0.  ,  78.7 ]])

Rotating tensors
=========

Rotating a Tensor in ZXZ Bunge Euler Convention
===============================================

To rotate a tensor using the ZXZ Bunge Euler convention, you perform a series of three rotations about specific axes. In this convention, the Euler angles :math:`(\alpha, \beta, \gamma)` represent rotations as follows:

1. **Rotation by :math:`\alpha`** (first angle) around the **Z-axis**.
2. **Rotation by :math:`\beta`** (second angle) around the **X-axis**.
3. **Rotation by :math:`\gamma`** (third angle) around the **Z-axis** (again).

The full transformation matrix for the ZXZ convention is the product of three individual rotation matrices:

Rotation Matrices
-----------------

**Rotation around the Z-axis** by angle :math:`\theta`:

.. math::

   R_Z(\theta) = \begin{pmatrix}
   \cos(\theta) & -\sin(\theta) & 0 \\
   \sin(\theta) & \cos(\theta) & 0 \\
   0 & 0 & 1
   \end{pmatrix}

**Rotation around the X-axis** by angle :math:`\phi`:

.. math::

   R_X(\phi) = \begin{pmatrix}
   1 & 0 & 0 \\
   0 & \cos(\phi) & -\sin(\phi) \\
   0 & \sin(\phi) & \cos(\phi)
   \end{pmatrix}

Full Rotation Matrix
--------------------

The full rotation matrix for the ZXZ convention is:

.. math::

   R = R_Z(\alpha) \cdot R_X(\beta) \cdot R_Z(\gamma)

Substitute the rotation matrices:

.. math::

   R = \begin{pmatrix}
   \cos(\alpha) & -\sin(\alpha) & 0 \\
   \sin(\alpha) & \cos(\alpha) & 0 \\
   0 & 0 & 1
   \end{pmatrix}
   \cdot
   \begin{pmatrix}
   1 & 0 & 0 \\
   0 & \cos(\beta) & -\sin(\beta) \\
   0 & \sin(\beta) & \cos(\beta)
   \end{pmatrix}
   \cdot
   \begin{pmatrix}
   \cos(\gamma) & -\sin(\gamma) & 0 \\
   \sin(\gamma) & \cos(\gamma) & 0 \\
   0 & 0 & 1
   \end{pmatrix}

Rotating a Tensor
-----------------

Given a second-order tensor :math:`T`, after applying the rotation matrix :math:`R`, the rotated tensor :math:`T'` can be obtained by:

.. math::

   T' = R \cdot T \cdot R^T

where :math:`R^T` is the transpose of the rotation matrix :math:`R`.

This process rotates the tensor according to the ZXZ Euler angles.

To rotate a tensor within santex, we can define alpha, beta and gamma in degrees, and then call the ``rotate_tensor()`` method within Tensor class as follows:

.. code-block:: python
    alpha = 10
    beta = 20
    gamma = 30

    rotated_forsterite = tensor.rotate_tensor(cijkl_forsterite, alpha, beta, gamma)
    voigt_rotated_forsterite = tensor.tensor_to_voigt(rotated_forsterite)

Now the same previous voigt tensor shown, is rotated, the new voigt notation of the rotated forsterite looks like:

.. code-block:: bash

    array([[254.40669486,  83.60527452,  74.25170921,   2.25426383,
          9.1917399 ,  32.19107185],
       [ 83.60527452, 229.09177813,  77.9967799 ,   3.16662267,
          4.15506621,  24.85704921],
       [ 74.25170921,  77.9967799 , 228.39399978,  -5.76818554,
          4.58001959,  -1.68082112],
       [  2.25426383,   3.16662267,  -5.76818554,  72.58824371,
          7.02836206,   6.26816834],
       [  9.1917399 ,   4.15506621,   4.58001959,   7.02836206,
         73.39492984,   5.79247342],
       [ 32.19107185,  24.85704921,  -1.68082112,   6.26816834,
          5.79247342,  93.02059007]])

Material analysis
----------

