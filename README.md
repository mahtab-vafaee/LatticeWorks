# LatticeWorks

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE) [![DOI](https://zenodo.org/badge/486929347.svg)](https://doi.org/10.5281/zenodo.13862475) 

This is a MATLAB toolbox for creating and customizing lattice structures, multi-morphology lattices, and lattices in different coordinates and arrangements. 

[![graphAbstract](docs/html/graphAbstract.png)](https://github.com/mahtab-vafaee/LatticeWorks/tree/main)

# Getting started
## Installation
To install this toolbox, simply run `install_me.m` found in the main project folder. 

Alternatively one can manually install by adding `lib` and `lib_ext` folders to the path. 

### Install dependancies
* Install the [GIBBON MATLAB toolbox](https://www.gibboncode.org/)
* Required MATLAB toolboxes:
	- Image Processing Toolbox
* ABAQUS, to run finite element simulations through ABAQUS
* [FEBio](https://www.febio.org/), to run finite element simulations through FEBio (see also configuration information with the GIBBON toolbox)
* [export_fig](https://github.com/altmany/export_fig), to fascilitate the creation of publication quality figures. 

## Running examples
Examples are contained in the `docs` folder.

# Citing LatticeWorks
If you make use of LatticeWorks, please cite the Materials and Design Journal paper: [![DOI](https://img.shields.io/badge/Materials%20and%20Design-10.1016/j.matdes.2024.113564-blue.svg)](https://doi.org/10.1016/j.matdes.2024.113564)
    
Vafaeefar, et al. (2025). LatticeWorks: An open-source MATLAB toolbox for nonuniform, gradient and multi-morphology lattice generation, and analysis. Materials & Design, 250, 113564, [https://doi.org/10.21105/joss.00506](https://doi.org/10.1016/j.matdes.2024.113564)

Please also cite Zenodo DOI [![DOI](https://zenodo.org/badge/486929347.svg)](https://doi.org/10.5281/zenodo.13862475) as a software citation.

# Applications
* Lattice generation for tissue engineering and scaffolds, biomedical devices, energy absorption, etc.
* Generating ready-to-print STL files for 3D printing of lattices. 
* Finite element analysis (FEA) using ABAQUS directly in the toolbox, as well as post-processing the results. `DEMO_0014_FEA_ABAQUS_Twisted_Cylindrical_Gyroid` is an example of FEA on a generated lattice structure, through ABAQUS directly in the toolbox.
* Mapping optimised nonuniform gradient lattices on distributed structural and mechanical properties, e.g. stiffness. `DEMO_0013_Mapping_Density_Distribution` maps nonuniform gradient gyroid on a density distribution field.
* Creat infill lattice structures within a closed surface, using different lattice types. `DEMO_0012_infill_STL_Lattice` is an example of this application on a vertebrae model.

# Contributing
We welcome submissions. If you'd like to contribute please file a pull request or post an issue. Thanks! 

# License <a name="License"></a>
LatticeWorks is provided under the [MIT license](https://github.com/mahtab-vafaee/LatticeWorks/blob/main/LICENSE).
