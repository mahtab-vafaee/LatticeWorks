---
title: 'LatticeWorks'
tags:
  - Matlab
  - lattice
  - TPMS
  - functionally graded lattice
  - multi-morphology
authors:
  - name: Mahtab Vafaeefar
    orcid: 0000-0002-0475-8586
    corresponding: true # (This is how to denote the corresponding author)
    affiliation: "1" # (Multiple affiliations must be quoted)
  - name: Kevin M. Moerman
    orcid: 0000-0003-3768-4269
    affiliation: "1, 2"
  - name: Ted J Vaughan
    orcid: 0000-0002-2219-7620
    affiliation: 1
affiliations:
 - name: Biomechanics Research Centre (BMEC), Schoole of Engineering, College of Science and Engineerig, University of Galway, Ireland
   index: 1
 - name: Griffifth University, Gold Coast, Australia
   index: 2
date: 28 August 2024
bibliography: paper.bib
---

# Summary

[LatticeWorks](https://github.com/mahtab-vafaee/LatticeWorks) is a MMATLAB toolbox for creating and customizing lattice structures, multi-morphology lattices, and lattices in different coordinates and arrangements. 

![A Graphical summary of the LatticeWorks toolbox](graphAbstract.png)

# Main Features and Usage

A summary of diffirent lattice configurations with LatticeWorks is presented here: 

* **Single-morphology lattices** which includes the introduction and generation of basic unit cells of TPMS, and then linearly graded lattice structures. 
* **Multi-morphology lattices**, describing different transition functions for generating a  gradient cell-type lattice structures (multi-morphology or hybrid).   
* **Boundary shapes**, to apply different boundary shapes on generated lattices (like cylindrical or spherical) rather than cubic samples, also to create arbitrary transition boundary shapes in multi-morphology lattices.  
* **Cell arrangement methods**, allows setting cells in different arrangements, rather than in Cartesian systems, and applying deformation tensors on the unit cells.  
* **Infill lattice structures**, featuring methods to volume infill any arbitrary surface with a lattice structure. 
* **Mapping density (non-)uniform graded lattices**, a general technique to map a tailored graded lattice on any density distribution field. This technique has been developed for covering the porous region with a solid shell as well. 
* **Finite Element Analysis**, providing an integrated FE tool to analyse designed lattices in terms of their mechanical properties and compare different lattice shapes. 
* **Designs for additive manufacturing prototypes**, preparing designed samples of functionally graded lattices for high resolution 3D printing.  

# Statement of need

Rapid advancement of additive manufacturing technology has expanded the design envelope for sophisticated lightweight lattice structures `[@Dong2022]`. Due to their highly tunable and multifunctional nature of functionally graded lattice structures, they are used extensively in different applications from biomedical to aerospace engineering `[@Perez2022; @Veloso2022]`. In biomedical applications these structures are used as bone implants and scaffolds `[@Vafaeefar2023; @Naghavi2023, @Zadpoor2019]`. In aerospace Lattice structures are frequently utilized in the aerospace industry along with topology optimization to provide lightweight designs `[@Veloso2022]`. Lattice structures are also used for energy absorption applications `[@Vafaeefar2024]`. Functionally graded lattice structures, either by unit cell porosity, cell size control, multi-morphology lattice structures, or even a combination of these techniques can be used to address the engineering specifications required `[@Al_ketan2020]`. 

Nevertheless, the literature clearly lacks a comprehensive set of tools to address different approaches for generating functionally graded lattice structures, more genrally, lattice structure in different configurations. Given the current significance and potential of functionally graded lattice structures for many applications, identifying techniques and tools that can control lattice geometries to construct functional gradients, and multi-morphology lattices is the main objective of this work. LatticeWorks offers the necessary tools to the researchers to easily customise a lattice structure, for a specific application.

# Acknowledgements

This project has received funding from the European Research Council (ERC) under the EU’s Horizon 2020 research and innovation program (Grant agreement No. 804108). This project has also received funding from the European Union’s Horizon Europe research and innovation programme under grant agreement No 101047008 (BIOMET4D). 

# References
        