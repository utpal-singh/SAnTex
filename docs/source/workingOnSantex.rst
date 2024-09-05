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


