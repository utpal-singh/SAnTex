---
#title: 'pide: an Earth material library for petrophysical calculations and synthetic modelling of Earth systems.'
title: 'SAGE: Seismic Anisotropy and Geodynamics - A Python-Based Tool for Comprehensive Seismic Anisotropy Calculation'
tags:
  - Python3
  - seismic anisotropy
  - petrophysics
  - geodynamics
  - seismic velocity
  - experimental petrology
  - geodynamics
  - geophysics

authors:
  - name: Utpal Singh
    orcid: 0000-0001-8304-5615
    equal-contrib: true
    corresponding: true 
    affiliation: 1 
  - name: Sinan Ozaydin
    orcid: 0000-0002-4532-9980
    affiliation: 1
  - name: Vasileios Chatzaras
    affiliation: 1
  - name: Patrice F. Rey
    affiliation: 1
  - name: Maria Constanza Manassero
    affiliation: 3

affiliations:
 - name: School of Geosciences, University of Sydney, Sydney, Australia.
   index: 1

date: 30 January 2024
bibliography: paper.bib
---

# Abstract
Seismic anisotropy, how seismic wave velocities vary in different directions, is a key parameter for understanding the structure and dynamics of the Earth's interior. However, measuring seismic anisotropy is challenging and requires accurate and reliable methods and tools. In this paper, we present SAGE: \textit{Seismic Anisotropy and Geodynamics}, a Python library that calculates seismic anisotropy from mineral phase, modal composition, and crystallographic preferred orientation data collected by electron backscatter diffraction (EBSD). SAGE incorporates mineral phase information, and a crystal stiffness tensor library accounting for the dependency of elasticity with pressure and temperature. The library also also allows for estimation of the volume of melt and its effect on seismic anisotropy. SAGE is an open-source and user-friendly python-based alternative to the MATLAB-based software MTEX, with some additional components with regards to geodynamics. We describe the design, implementation, and validation of SAGE against previous studies and theoretical models. SAGE is the first library for seismic anisotropy calculation that accounts for the effect of pressure and temperature on elastic stiffness tensors.SAGE can produce reliable and consistent results that agree with previous studies, scientific tools and theoretical models.

# Statement of need
Our understanding of the dependency between seismic anisotropy and microstructures relies on laboratory measurements and calculations of seismic properties based on rock microfabric. The microfabric properties such as, crystallographic preferred orientation, grain sizes and distribution of microcracks, caninfluence both seismic velocity and anisotropy  (Almqvist and Mainprice 2017).

Seismic anisotropy calculations based on EBSD data and standard reference crystal stiffness tensors, have proven indispensable in understanding seismic anisotropy (Citations). However, advancements in density functional theory and laboratory simulations under varying pressure and temperature conditions highlight significant deviations in crystal stiffness tensors, which are underutilized (Kumazawa and Andersen 1969, Kumazawa 1969, Walker 2012, Qian et al. 2017, Su et al. 2021). Given the range and diversity of geothermal settings within the Earth, it becomes imperative to acknowledge and incorporate these variations. Consequently, there is a need for a seismic analysis tool that accommodates the pressure and temperature responses of crystal stiffness tensors, ensuring a more accurate portrayal of seismic wave velocity and anisotropy.

Moreover, rocks may undergo partial melting, which impacts strongly on their elastic properties, influencing both seismic wave velocity and anisotropy (e.g., Lee et al 2017). Hence, we need a tool that not only considers the impact of pressure and temperature on crystal properties but also incorporates the effects of melt presence on seismic properties. 

Adding another layer of complexity, variations in phase density under changing pressure and temperature conditions further impact seismic wave velocity. We have implemented a Python module that calculates these changes in density and seismic wave speeds if the phases are considered isotropic (Hacker et al 2003).

While existing solutions predominantly utilize MATLAB for such calculations, its proprietary status present barriers to widespread accessibility. Consequently, and to facilitate seamless and open-source access to seismic speed and anisotropy analyses, we propose the development of a dedicated Python module. This module, named SAGE (Seismic Anisotropy and Geodynamics), serves as a comprehensive, freely accessible package for researchers, incorporating functionalities for pressure and temperature, melt, and density calculations. By making these essential tools available in Python, we aim to simplify seismic analysis, and accessible platform for the broader scientific community.

# Package Summary

`SAGE` (Seismic Anisotropy and Geodynamics) is a comprehensive Python-based package designed for calculation of seismic anisotropy by providing a versatile, open-source toolkit for researchers. SAGE caters to the dynamic nature of crystal properties, integrating functionalities for seismic anisotropy calculations, pressure-temperature effects, melt presence, and phase density variations.

Within `sage`, `EBSD` is the class which deals with ebsd based calculations. The relevant functions for users are `EBSD.load(filename)`, `EBSD.filterByPhase(phases)`, `EBSD.orientations()`, `EBSD.plot()`

`Tensor` is the class which deals with tensor calculations. `Tensor.rotate()` will rotate the orientations obtained from `EBSD.orientations()` to sample reference frame. `Tensor.voigtToTensor()` converts a voigt 6*6 tensor to 3*3*3*3 tensor notation, and `Tensor.tensorToVoigt()` converts a 3*3*3*3 tensor notation to voigt notation.

`Material` is the class that outputs an user defined material object. The user can input if they want isotropic velocities, or anisotropic velocities, they get the library of predefined values of thermodynamic constants, stiffness tensors. The user also gets to decide the pressure, temperature and melt percentage if they want to include that in their material class. The user can choose their own density or get the density calculated from `Material.thermodynamics.calculateDensity()` 

`Anisotropy` is the class that calculates seismic anisotropy. The `Tensor.rotate()` object and `Material` object is taken as input into `Anisotropy.getVelocity()` and `Anisotropy.getAnisotropyPlot()`

`Isotropy` is the class that calculates vp, vs, density, bulk modulus from library values at any chosen pressure and temperature.

![Workflow Chart for SAGE \label{fig:SAGE_wflow}](docs/images/SAGE_workflow.png)

# Acknowledgements
This study is supported by the Australian Research Council (ARC)Linkage Grant #Grantnumber and ARC DP Grant #Grantnumber. 

# References
