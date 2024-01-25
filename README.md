# OpenBoltRF
## An open-source workflow to analyse load transfer for geometrically imperfect bolted ring-flanges in offshore wind turbines.

OpenBoltRF is an open-source workflow to analyse the load transfer for geometrically imperfect bolted ring-flanges in offshore wind turbine support structures. It is intended as a lightweight, portable package for users familiar with FreeCAD, Salome Meca / Code Aster, Paraview, and Python. The package is introduced in a paper from the TORQUE 2024 conference (<span style="color:blue">some *blue* text</span>.) 

[![DOI](https://zenodo.org/badge/683682811.svg)](https://zenodo.org/badge/latestdoi/683682811)

# Structure
The package contains a script for controlling all actions, _main.py_. 

Within _main.py_, a datafile specifying the parameters of the flange (dimensions, material properties), _Flange_params.py_, is imported. A setup file specifying the particulars of the project (filenames, paths etc.), _setup.py_, is also imported by _main.py_. The values for variables in both _Flange_params.py_ and _setup.py_ may be configured by the user, but the variables names must not be altered.

With the flange parameters and project details specified the execution of the workflow may be completed. The general order is as follows:-

  - The geometry is generated using the _geom_ module. _geom_ calls **FreeCAD**, and executes a script in the **FreeCAD Python console** to generate .stp files for the flange under study.
 
  - A mesh is generated from the .stp files, using the _mesh_ module. _mesh_ calls **Salome Meca**, and executes a script in the **Salome Meca Python console** to generate the meshes.
 
  - A simulation is executed with the generated mesh, using the _sim_ module. _sim_ calls **Code Aster**, and executes a pre-defined .comm file to run the simulation.
    - The .comm file specifies the particular of the **Code Aster** FEM simulation, and is commented accordingly. A pre-existing knowledge of Code Aster is required.
    - The simulation is executed using an _export_ file, which specifies simulation details such as number of processor cores, simulation memory limit, scratch location, output file location etc. Such parameters are configured in the _setup_ file dictated above.
    - Proceeding past this point in the workflow is contingent on the succesful convergence of the _Code Aster_ simulation.
    
  - Post processing is executed on a successful=ly converged FEM simulation result, using the _post_pro_ module. _post_pro_ calls **Salome Meca**, and executes a script in the **Paravis Python console** to generate summary results for the simulation.
  
  - Further analysis and visualisation TBA. 

# Pre-requisites
The use of OpenBoltRF assumes execution on Ubuntu (tested on 18.04.6 LTS, no guarantee on other versions), with access to the following pre-requisites:- 
  - FreeCAD 0.21.0     or later
  - Salome Meca 2021   or later
    - Singularity container installation, including Paravis
  - Python 3.8.5       or later
    - matplotlib>=3.3.2
    - numpy>=1.19.2
    - pandas>=1.2.1
   
# Instructions for use
At this stage OpenBoltRF is not designed as a Python package. It is only designed to be downloaded to a relevant project folder, and be called directly in that project folder.

  - Download the repository to '/path/to/your/project'
    
  - Edit the _setup.py_ file to tailor to your project path and simulation requirements
    
  - Edit the _Flange_params.py_ datafile to tailor to your bolted ring-flange
    
  - Run the _main.py_ script to generate a result

# Contact
Please email jack.jorgensen@research.uwa.edu.au for any queries related to the workflow. Any queries specific to the tools used (FreeCAD, Salome Meca / Code Aster, & Paraview) should be directed toward the relevant software communities.
