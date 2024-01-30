#!/usr/bin/env python
# -*- coding: utf-8 -*-

# region postp_pv.py 
## Description 
## OpenBoltRF postp_pv.py file, psot processing FEM simulation results to extract relevant stress state data in bolted ring-flanges using Paraview
__author__ = "Jack Jorgensen"
__license__ = "GNU GPLv3"
__version__ = "Per Git commit"
__maintainer__ = "Jack Jorgensen"
__email__ = "jack.jorgensen@research.uwa.edu.au"
__status__ = "Development"
# endregion

# region Import packages
    # Paraview packages to use Paravis inside Salome Meca container
from paraview.simple import *

    # System packages
import os
import subprocess
import platform
from pathlib import Path
import sys
import numpy as np
import math
import argparse
import csv

# endregion

# region Set paths

parser = argparse.ArgumentParser()

parser.add_argument('-out', action='store', dest='out')
parser.add_argument('-run', action='store', dest='run')
parser.add_argument('-path_pre_folder', action='store', dest='path_pre_folder')

args = parser.parse_args()
out_name = args.out
run_path = args.run
path_pre_folder = args.path_pre_folder

sys.path.append(path_pre_folder)

from flange_params import * # Import flange parameters to allow calc set-up

# endregion

# region Post Process

outrmed = MEDReader(registrationName=out_name, FileName=path_pre_folder + run_path + out_name)
outrmed.AllArrays = ['TS0/00000001/ComSup1/Load_steSIEQ_NOEU@@][@@P1']
outrmed.GenerateVectors = 1
timesteps = outrmed.TimestepValues

# endregion

# region Create probe locations
probeLocations = {}
probetime = {}

for i in range(0, len(Probes)):
	probeLocations[i] = ProbeLocation(registrationName='ProbeLocation'+str(i), Input=outrmed, ProbeType='Fixed Radius Point Source')
	probeLocations[i].ProbeType.Center = Probes[i]
	
	
# endregion

# region Get data over quasi-static simulation
groupTimeSteps = {}

for i in range(0, len(Probes)):
	groupTimeSteps[i] = GroupTimeSteps(registrationName='GroupTimeSteps1'+str(i), Input=probeLocations[i])
# endregion

# region Save .csv of probe data
	# Save data
for i in range(0, len(Probes)):
    SaveData(path_pre_folder + run_path+'data_'+out_name+'_'+str(i)+'.csv', proxy=groupTimeSteps[i], ChooseArraysToWrite=1,
        PointDataArrays=['Load_steSIEQ_NOEU', 'Load_steSIEQ_NOEU_Vector'],
        AddMetaData=1,
        AddTime=1)

	# Append time
with open(path_pre_folder + run_path+'time_'+out_name+'.csv', 'w', newline='') as csvfile:
    timewriter = csv.writer(csvfile, delimiter=',',
                            quotechar='|', quoting=csv.QUOTE_MINIMAL)
    timewriter.writerow(timesteps)

# endregion