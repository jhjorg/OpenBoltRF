#!/usr/bin/env python
# -*- coding: utf-8 -*-

# region mesh_imperf.py 
## Description 
## OpenBoltRF mesh_imperf.py file, adjusting .med mesh files for imperfections
__author__ = "Jack Jorgensen"
__license__ = "GNU GPLv3"
__version__ = "Per Git commit"
__maintainer__ = "Jack Jorgensen"
__email__ = "jack.jorgensen@research.uwa.edu.au"
__status__ = "Development"
# endregion

# region Import packages
    # Salome packages
import salome_notebook
import  SMESH, SALOMEDS
from salome.smesh import smeshBuilder

    # System packages
import sys
import numpy as np
import math
import platform

# endregion

 # region Set path
path_pre_folder = sys.argv[4]
# endregion

# region Flange Parameters & project definition

sys.path.append(path_pre_folder)

from flange_params import * # Import flange parameters to allow calc set-up

project_Number = sys.argv[2]
job_name = sys.argv[3]
run_num = sys.argv[7]

file_Name = "Flange_" + project_Number
run_filepath = path_pre_folder + "Runs/" + job_name

# endregion

# region Node moves

    # Load mesh
smesh = smeshBuilder.New()

([Compound_Mesh_1], status) = smesh.CreateMeshesFromMED(run_filepath + '/Compound_Mesh_'+file_Name+'.med')

    # Imperfection angle & height
alpha = float(sys.argv[5])
k_height = float(sys.argv[6])

    # Get nodes - top & bottom flange
Group_names = Compound_Mesh_1.GetGroupNames()

for i in range(0, len(Group_names)):
    if Group_names[i] == 'T_Flange_con':
       index_TF = i
       node_it_TF = Compound_Mesh_1.GetGroups()[i].GetListOfID()
    elif Group_names[i] == 'B_Flange_con':
       index_TF = i
       node_it_BF = Compound_Mesh_1.GetGroups()[i].GetListOfID()

    # Move nodes
tolerance = 1e-6 # Set tolerance for nodes
selected_nodes = [] # Initilise empty list
radius = [] # List for radii
theta = [] # List for angle

    # Loop to check if perfect flange
if alpha != 0:
            # Top flange
    for node in node_it_TF:
        node_coords = np.array(Compound_Mesh_1.GetNodeXYZ(node))
        if (abs(abs(node_coords[2]) - (flange_height+flange_lip+shell_height)/1000) <= tolerance):
            selected_nodes.append(node)
            radius.append((node_coords[0]**2+node_coords[1]**2)**0.5)
            theta.append((180/np.pi)*np.arctan2(node_coords[1], node_coords[0]))
            theta[-1] = abs(theta[-1]-180)        
            if theta[-1] > (alpha/2):
                amp_up_node = 0
            else:
                amp_up_node = ((k_height/4)*np.sin(2*np.pi*(1/alpha)*theta[-1]+np.pi/2)+(k_height/4))/1000
            isDone = Compound_Mesh_1.MoveNode( node, node_coords[0], node_coords[1], node_coords[2]+amp_up_node )
            # Bottom flange
    for node in node_it_BF:
        node_coords = np.array(Compound_Mesh_1.GetNodeXYZ(node))
        if (abs(abs(node_coords[2]) - (flange_height+flange_lip+shell_height)/1000) <= tolerance):
            selected_nodes.append(node)
            radius.append((node_coords[0]**2+node_coords[1]**2)**0.5)
            theta.append((180/np.pi)*np.arctan2(node_coords[1], node_coords[0]))
            theta[-1] = abs(theta[-1]-180)      
            if theta[-1] > (alpha/2):
                amp_up_node = 0
            else:
                amp_up_node = ((k_height/4)*np.sin(2*np.pi*(1/alpha)*theta[-1]+np.pi/2)+(k_height/4))/1000
            isDone = Compound_Mesh_1.MoveNode( node, node_coords[0], node_coords[1], node_coords[2]-amp_up_node )

# endregion

# region Save mesh
try:
  Compound_Mesh_1.ExportMED(run_filepath + f"/Run_{str(run_num)}" + '/Compound_Mesh_'+file_Name+'_'+str(run_num)+'.med',auto_groups=0,minor=40,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')        

# endregion