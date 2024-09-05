.. santex documentation master file, created by
   sphinx-quickstart on Thu Sep  5 21:41:54 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to santex's documentation!
==================================

`Seismic Anisotropy from Texture` (SAnTex) is a Python library which calculates the full elastic tensor of rocks from modal mineral composition, crystallographic orientation, and a crystal stiffness tensor catalogue that accounts for the dependency of elasticity with pressure and temperature.


----------

Features
--------

* Pre-processing and cleaning of EBSD data: SAnTex facilitates the processing and cleaning of EBSD data. To enhance data completeness, SAnTex offers the option to fill not-indexed pixels or indexed pixels removed during the cleaning process, using machine learning techniques.
* Tensor operations: Tensor conversions between Voigt matrix and full stiffness tensors, as well as rotations based on euler angles.
* Material analysis: SAnTex provides a catalogue of minerals, users can load the catalogue and can either utilise them to load stiffness tensors for their EBSD phases or make a modal rock and can calculate seismic anisotropy for the modal rock.
* Seismic Anisotropy: SAnTex performs calculations of seismic anisotropy at a range of pressure and temperature conditions. It also offers visualisation capabilities, allowing users to view the calculated seismic anisotropy in 2D and 3D plots.
* Isotropic velocities: Calculates isotropic seismic wave velocities (Vp, Vs and vbulk), isothermal bulk modulus, and density at elevated temperatures and pressures (Hacker & Abers, 2004).


Installation
-------------

install.rst

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   santex.ebsd.rst
   santex.material.rst     
   santex.tensor.rst               
   santex.anisotropy.rst   
   santex.isotropy.rst     
   santex.rst
   modules

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
