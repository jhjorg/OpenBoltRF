<img src= 'https://github.com/jhjorg/OpenBoltRF/assets/70005890/f61b2e8c-bf70-48b4-bbbc-d7468b080c32' width = '400'>

# OpenBoltRF
## An open-source workflow to analyse load transfer for geometrically imperfect bolted ring-flanges in offshore wind turbines.

OpenBoltRF is an open-source workflow to analyse the load transfer for geometrically imperfect bolted ring-flanges in offshore wind turbine support structures. It is intended as a lightweight, portable package for users familiar with FreeCAD, Salome Meca / Code Aster, Paraview, and Python. The package is introduced in a paper from the TORQUE 2024 conference entitled _An open-source analysis workflow for geometrically imperfect bolted ring-flanges in wind turbine support structures_ - ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png) (LINK TO PAPER HERE) ![#f03c15](https://placehold.co/15x15/f03c15/f03c15.png)

[![DOI](https://zenodo.org/badge/683682811.svg)](https://zenodo.org/badge/latestdoi/683682811)

We use four open-source packages in our workflow:- the Python programming language, the FreeCAD 3D parametric modelling package, the Code Aster / Salome-Meca FEM simulation environment, and the ParaView post-processing package. 


# Structure

The workflow is summarised in the figure below. 

  - All actions are controlled by a main Python script, _main.py_.

  - A Python script with all flange dimensions is prepared, named _flange_params.py_.

  - A Python script, _geom_frc.py_, is executed via a subprocess command at the FreeCAD command line, to create a parametric flange geometry based on _flange_params.py_. The flange geometry is saved as 4 .stp files - the top & bottom flange + shells, the bolt, and the nut.

  - A Python script, _mesh_sm.py_, is executed via a subprocess command at the Salome Meca command line, to create a compound mesh object of the .stp geometries. The geometry imperfections are implemented in with the script _mesh_imperf.py_, and can be explicitly specified with height / angle, or can be specified as random according to the work of Buchholz & Seidel [1].
  
  - The Code Aster FEM solver is called via a subprocess command to run a simulation based on defined _export_ and _.comm_ files, and using the generated mesh. The _.comm_ file defines the simulation features, and the _export_ file defines the simulation management parameters (which _.comm_ file to use, number of processor cores, simulation memory limit, scratch location, output file location etc.).
  
  - A Python script, _postp_pv.py_, is executed via a subprocess command at the Paraview command line inside Salome Meca, to post-process FEM simulation results. Useful results, including the decomposed axial and bending loads in particular bolts on the flange are calculated.

  - A Python script, _fatigue.py_, is executed to calculate fatigue damage from a given load cycle, given the load transfer function determined via the FEM simulations.  

![Workflow](https://github.com/jhjorg/OpenBoltRF/assets/70005890/29766f6e-cfc9-46c1-9e97-650cdcbe18b3)

Mesh and parametric studies can be performed with the provided looping structure in the _main.py_ script. 

# Pre-requisites
The use of OpenBoltRF assumes execution on a Linux OS (tested on Ubuntu 18.04.6 LTS, no guarantee on other versions), with access to the following pre-requisites:- 
  - FreeCAD 0.21.0     or later (https://github.com/FreeCAD/FreeCAD-Bundle/releases/tag/0.21.0)
  - Salome Meca 2021   or later (https://code-aster.org/FICHIERS/singularity/salome_meca-lgpl-2021.0.0-0-20210601-scibian-9.sif)
    - Singularity container installation, including Paravis distribution of Paraview
  - Python 3.8.5       or later 

Instructions on installation of these tools can be found from respective websites. Note, not tested on MPI version of Code Aster / Salome Meca.

The Python packages required for execution on the _main.py_ and _fatigue.py_ scripts are collected in _requirements.txt_. The extra Python packages used in other scripts (for example _salome_, _GEOM_ etc. in _mesh_sm.py_) are application specific, and should be installed in the packaged Python distribution internal to the applications (FreeCAD, Salome Meca and ParaView).
   
# Instructions for use
At this stage OpenBoltRF is not designed as a Python package. It is only designed to be downloaded to a relevant project folder, and be called directly in that project folder.

  - Download the repository to '/path/to/your/project'
    
  - Edit the _main.py_ file to tailor to your project path
    
  - Edit the _flange_params.py_ datafile to tailor to your bolted ring-flange

  - Edit the _.comm_ and _export_ files to tailor to your simulation settings (no. of processor, RAM, scratch location etc.)
    
  - Run the _main.py_ script to generate a result

# Contact
Please email jack.jorgensen@research.uwa.edu.au for any queries related to the workflow. Any queries specific to the tools used (FreeCAD, Salome Meca / Code Aster, & Paraview) should be directed toward the relevant software communities.

# References
[1] Buchholz A and Seidel M 2023, "Gap height prediction for bolted ring flange connections based on measurements", _Steel Construction_ 16 114–126
