#!/usr/bin/env python
# -*- coding: utf-8 -*-

# region main.py 
## Description 
## OpenBoltRF main file, calling other modules
__author__ = "Jack Jorgensen"
__license__ = "GNU GPLv3"
__version__ = "Per Git commit"
__maintainer__ = "Jack Jorgensen"
__email__ = "jack.jorgensen@research.uwa.edu.au"
__status__ = "Development"
# endregion

# region Import required packages and modules
    # Packages
import numpy as np
import math
from matplotlib import pyplot as plt
import os
import subprocess
import platform
from pathlib import Path
import sys
import pandas as pd
import scienceplots
from scipy.stats import lognorm

    # Set path
path_pre = "/path/to/your/project/"
                
folder = "/OpenBoltRF/"
path_pre_folder = r"{}".format(path_pre+folder)

    # Flange Parameters
sys.path.append(path_pre_folder)

from flange_params import * # Import flange parameters to allow calc set-up

# endregion

# region Set paths for running

    # Set the path to your Singularity container
container_path = "/set/path/to/SalomeMeca/container.sif"

    # Set the job type - "Mesh" for mesh sweep, "Sweep" for parametric sweep
job_type = "Sweep"

    # Set the job name - describe what you doing the runs for
job_name = "Moment_runs"

    # Set paths to run files
run_filepath = path_pre_folder + "Runs/" + job_name
if not os.path.exists(run_filepath):
    os.makedirs(run_filepath)

# endregion

# region Create geometry
project_Number = "Project_Number"
file_Name = "Flange_" + project_Number

freeCAD_appImage = "/path/to/FreeCAD/.AppImage"

geometry_command = freeCAD_appImage + " freecadcmd " + path_pre_folder + "geom_frc.py --pass '" +project_Number+"' '"+job_name+"' '"+path_pre_folder+"'"

if not os.path.isfile(run_filepath+"/" + file_Name + "_Top_Flange_dual_mat_multi.step"):
    subprocess.call(geometry_command, timeout=600, shell=True)
# endregion

# region Create mesh without imperfections
    # Folder mapping - # May be necessary to map working folder to Salome Meca singularity container.
mapping  = 1

if mapping == 1:
    mapped_folder = "/mnt"  
    path_pre_mapped = r"{}".format(mapped_folder+folder)
    map_command = " -B " + path_pre + ":" + mapped_folder + " " 
elif mapping == 0:
    mapped_folder = path_pre  
    path_pre_mapped = r"{}".format(mapped_folder+folder)
    map_command = " " 

    # Mesh command
mesh_command = "singularity run" + map_command + container_path+" shell -- python "+path_pre_mapped+ "mesh_sm.py --pass '" +project_Number+"' '"+job_name+"' '"+path_pre_mapped+"'"

    # Run Salome Meca container, and execute mesh_command at Python command line inside container
if not os.path.isfile(run_filepath + '/Compound_Mesh_'+file_Name+'.med'):
    subprocess.call(mesh_command, timeout=600, shell=True)
# endregion
    
# region Implement mesh imperfections
    # Imperfection definitions
imperf_def = "Random"
if imperf_def == "Random":
    if not os.path.isfile(run_filepath + '/Imperfs_'+file_Name+'.txt'):
            # Run random sampling - per Buchholz & Seidel 2023
        Angle = 30      # Set angle of imperfection
        num_samps = 10  # Set number of samples
        l_k = (Angle/360)*(np.pi*flange_out_dia/1000)
        E_k = 0.01*(l_k**2) + 0.2*l_k
        COV_k = (l_k)**(-1.8)+0.4
        sig_k = (np.log(COV_k**2+1))**0.5
        mu_k = np.log(E_k)-0.5*sig_k**2
        dist=lognorm(s=sig_k,scale=math.exp(mu_k))
        Imperf_height = dist.rvs(size=num_samps)
        Imperf_angle = [Angle] * num_samps
        np.savetxt(run_filepath + '/Imperfs_'+file_Name+'.txt', [Imperf_height, Imperf_angle])
    else:
        imperfs = np.loadtxt(run_filepath + '/Imperfs_'+file_Name+'.txt')
        Imperf_height = list(imperfs[0,:])
        Imperf_angle = list(imperfs[1,:])

else:
    if not os.path.isfile(run_filepath + '/Imperfs_'+file_Name+'.txt'):
        Imperf_angle  = [0, 30, 120, 25]    # Set defined imperfection angle
        Imperf_height = [0, 2, 2, 2]        # Set defined imperfection height
        np.savetxt(run_filepath + '/Imperfs_'+file_Name+'.txt', [Imperf_height, Imperf_angle])
    else:
        imperfs = np.loadtxt(run_filepath + '/Imperfs_'+file_Name+'.txt')
        Imperf_height = list(imperfs[0,:])
        Imperf_angle = list(imperfs[1,:])

    # Implement imperfection on new meshes
for i in range(0,len(Imperf_angle)):
        # Set-up folder locations for run files
    newpath_host = path_pre_folder + "Runs/" + job_name + f"/Run_{str(i)}"
    newpath_mapped = path_pre_mapped + "Runs/" + job_name + f"/Run_{str(i)}"

        # Create new run folders (if needed)
    if not os.path.exists(newpath_host):
        os.makedirs(newpath_host)

        # Imperfection command
    imperf_command = "singularity run" + map_command + container_path+" shell -- python "+path_pre_mapped+ "mesh_imperf.py --pass '" +project_Number+"' '"+job_name+"' '"+path_pre_mapped+"' '"+str(Imperf_angle[i])+"' '"+str(Imperf_height[i])+"' '"+str(i)+"'"

        # Run imperfection command
    if not os.path.isfile(newpath_host + '/Compound_Mesh_'+file_Name+'_'+str(i)+'.med'):
        subprocess.call(imperf_command, timeout=600, shell=True)  

# endregion

# region Code Aster loop
    
    # Set the path to your export file
export_file_path = path_pre_folder + "export"

    # Set the path to your comm file
comm_file_path = path_pre_folder + "Flange_run_multi.comm"

    # Set tmp file location for new runs
tmp_folder_path = mapped_folder + "/tmp"

    # Set path to meshes
mesh_file_paths = []
for i in range(0,len(Imperf_angle)):
    mesh_file_paths.append(path_pre_mapped + 'Runs/' + job_name + f"/Run_{str(i)}" + '/Compound_Mesh_'+file_Name+'_'+str(i)+'.med')

    # Find lines in export file
with open(export_file_path, "r") as f:
    export_file_lines = f.readlines()

        # Find the line with the .med file path, studyid, and .rmed results files
mesh_file_line = [line for line in export_file_lines if ".med" in line][0]
result_file_line = [line for line in export_file_lines if ".rmed" in line]
studyid_line = [line for line in export_file_lines if "studyid" in line][0]
comm_file_line = [line for line in export_file_lines if "comm" in line][0]
mess_file_line = [line for line in export_file_lines if "mess" in line][0]
nomjob_line = [line for line in export_file_lines if "nomjob" in line][0]
base_line = [line for line in export_file_lines if "R base" in line][0]

    # Find lines in comm file
with open(comm_file_path, "r") as f:
    comm_file_lines = f.readlines()

        # Find the line with the .med file path, studyid, and .rmed results files
Dim_ind = ["sig_F_v_b", "sig_F_v_n", "sig_F_z", "M_y"]  
preload_b_line = [line for line in comm_file_lines if Dim_ind[0] in line][0]
preload_n_line = [line for line in comm_file_lines if Dim_ind[1] in line][0]
if load_switch == 0:
    load_line = [line for line in comm_file_lines if Dim_ind[2] in line][0]
else:
    load_line = [line for line in comm_file_lines if Dim_ind[3] in line][0]

    # Create a set of new export files (and, optinally, .comm files) that refer to the different runs in the folder, delete studyid line, and point to correct tmp folder location
        # Set-up number of iterations
if job_type == "Single":
    job_num=1
else:
    job_num = len(mesh_file_paths)

for i in range(0, job_num):

    # Set mesh file paths - iterate if iteration, static if other
    if job_type == "Mesh":
        mesh_file_path = mesh_file_paths[i]
    elif job_type != "Mesh":
        mesh_file_path = mesh_file_paths[0]

    # Set tmp folder path for each run
    tmp_run_path = tmp_folder_path + '/Run_' + str(i)

    # Set-up folder locations for run files
    newpath_host = path_pre_folder + "Runs/" + job_name + f"/Run_{str(i)}"
    newpath_mapped = path_pre_mapped + "Runs/" + job_name + f"/Run_{str(i)}"

    if not os.path.exists(newpath_host):
        os.makedirs(newpath_host)

        # Read & write new .comm file
    preload_bolt = str(-1 * F_v / (math.pi * bolt_stress_dia/1000 * (nut_height+washer_height)/1000))
    preload_nut = str(1 * F_v / (math.pi * bolt_stress_dia/1000 * (nut_height+washer_height)/1000))
    if load_switch == 0:    
        shell_load = str(sig_F_z)
    else:
        shell_load = str(M_y)

    with open(newpath_host + '/Flange_run_multi_'+str(i)+'.comm','w') as outfile:
        with open(comm_file_path, "r") as f:
            a=f.read()
            if preload_b_line in a:
                    b=a.replace(preload_b_line, '    preload[i] = AFFE_CHAR_MECA(FORCE_FACE=(_F(FZ='+preload_bolt+',\n')
            if preload_n_line in b:
                    c=b.replace(preload_n_line, '                                        _F(FZ='+preload_nut+',\n')
            if load_line in c:
                    if load_switch == 0: 
                        d=c.replace(load_line, '    Seg_Load = AFFE_CHAR_MECA(FORCE_FACE=_F(FZ='+shell_load+',\n')
                    else:
                        d=c.replace(load_line, "                           FORCE_NODALE=_F(GROUP_NO='Moment_node',MY="+shell_load+",),);\n")

            outfile.write(d)

        # Read in the export file
    with open(newpath_host + '/export_'+str(i),'w') as outfile:
        with open(export_file_path, "r") as f:
            a=f.read()
            if mesh_file_line in a:
                #print (a)
                b=a.replace(mesh_file_line,'F libr ' + mesh_file_path + ' D  4\n')
                
            if result_file_line[0] in b:
                #print (a)
                c = b.replace(result_file_line[0],'F rmed ' + newpath_mapped + '/Out1_' + str(i) + '.rmed R  3\n')
                
            if base_line in c:
                #print (a)
                d = c.replace(base_line,'\n')
                
            if studyid_line in d:
                e = d.replace(studyid_line,'')

            if comm_file_line in e:
                f = e.replace(comm_file_line, 'F comm ' + newpath_mapped + '/Flange_run_multi_'+str(i)+'.comm' + ' D  1\n')
            
            if mess_file_line in f:
                g = f.replace(mess_file_line, 'F mess ' + newpath_mapped + '/message R  6\n')
 
            if nomjob_line in g:
                h = g.replace(nomjob_line, 'P nomjob RunCase_' + mesh_file_path + '\n')            
            outfile.write(h)
    
    with open(newpath_host + '/export_'+str(i), 'a') as file:
        file.write('P proxy_dir ' + tmp_folder_path +'\n\nP rep_trav ' + tmp_run_path + '\n')

    # Loop over jobs, and excute Code Aster runs for each
for i in range(0, job_num):

    print(i)
    # Reset run path
    newpath_mapped = path_pre_folder + "Runs/" + job_name + f"/Run_{str(i)}"

    # Set the command you want to run inside the container using as_run with the --run option
    command = f"-- as_run --run " + newpath_mapped + '/export_'+str(i)

    # Build the subprocess command to call the Singularity container
    subprocess_command = f"singularity run --pid{map_command}{container_path} shell {command}"

    # Call the subprocess command, and catch any errors 
    try:
            # Run command
        subprocess.call("exec "+subprocess_command, timeout=172800.0, shell=True)
            # Write error log
        error_log = "Succesful run - see output files"
            # Save error log
        with open(path_pre_folder + "Runs/" + job_name + f"/Run_{str(i)}/error_log.txt","w+") as f:
            f.write(error_log)
            # Catch timeout
    except subprocess.TimeoutExpired:
        error_log = "Run timeout expired - please re-try this iteration"
        print(error_log)
        with open(path_pre_folder + "Runs/" + job_name + f"/Run_{str(i)}/error_log.txt","w+") as f:
            f.write(error_log)

print("Iterations finished")

# endregion       

# region Paraview post-processing to extract the stresses required

for i in range(0, job_num):
    command = f"""shell -- python {path_pre_mapped}postp_pv.py -out 'Out1_{i}.rmed' -run 'Runs/{job_name}/Run_{i}/' -path_pre_folder '{path_pre_mapped}'"""

    # Build the subprocess command to call the Singularity container
    subprocess_command = f"""singularity run -B /data:/mnt {container_path} {command}"""

    # Call the subprocess command, and catch any errors 
    try:
            # Run command
        subprocess.call(subprocess_command, timeout=600, shell=True)
            # Write error log
        error_log = "Succesful run - see output files"
            # Save error log
        with open(path_pre_folder + "Runs/" + job_name + f"/Run_{str(i)}/error_log_PV.txt","w+") as f:
            f.write(error_log)
            # Catch timeout
    except subprocess.TimeoutExpired:
        error_log = "Run timeout expired - please re-try this iteration"
        with open(path_pre_folder + "Runs/" + job_name + f"/Run_{str(i)}/error_log_PV.txt","w+") as f:
            f.write(error_log)


# endregion

# region Comparison of results 

    # Import .csv files
probeData = {}

for j in range(0, job_num):
    probeData["Run_"+str(j)] = {}
    probeData["Run_"+str(j)]["Time"] = pd.read_csv(path_pre_folder + "Runs/" + job_name +"/Run_"+ str(j)+"/time_Out1_"+str(j)+".rmed.csv", header=None).T.round(decimals=2)

    for i in range(0, len(Probes)):
        probeData["Run_"+str(j)][i] = pd.read_csv(path_pre_folder + "Runs/" + job_name +"/Run_"+ str(j)+"/data_Out1_"+str(j)+".rmed_"+str(i)+".csv")
        probeData["Run_"+str(j)][i].index = (probeData["Run_"+str(j)]["Time"][0].tolist())
        probeData["Run_"+str(j)][i].index = pd.to_datetime(probeData["Run_"+str(j)][i].index, unit='s')

    # Check all probe points and scale if necesary
    if (probeData["Run_"+str(j)][0]['Load_steSIEQ_NOEU_Vector:0'][0] == 0):
        probeData["Run_"+str(j)][0]['Load_steSIEQ_NOEU_Vector:0'] = (probeData["Run_"+str(j)][1]['Load_steSIEQ_NOEU_Vector:0']) - (probeData["Run_"+str(j)][2]['Load_steSIEQ_NOEU_Vector:0'] - probeData["Run_"+str(j)][1]['Load_steSIEQ_NOEU_Vector:0'])

    # Calculate tension and bending components of bolt load at each bolt station
sig_tens = {}
sig_bend = {}
M_bolt = {}
F_bolt = {}
sig_tens_As = {}
sig_bend_As = {}
delta_sig_total = {}
delta_sig_tens = {}
delta_sig_bend = {}

for i in range(0, job_num):
    sig_tens[i] = (probeData["Run_"+str(i)][0]['Load_steSIEQ_NOEU_Vector:0'].resample('0.2S').interpolate() + probeData["Run_"+str(i)][len(probe_locs)-1]['Load_steSIEQ_NOEU_Vector:0'].resample('0.2S').interpolate())/2 # Pa, Bolt tension stress at shank
    sig_bend[i] = (probeData["Run_"+str(i)][len(probe_locs)-1]['Load_steSIEQ_NOEU_Vector:0'].resample('0.2S').interpolate() - probeData["Run_"+str(i)][0]['Load_steSIEQ_NOEU_Vector:0'].resample('0.2S').interpolate())/2 # Pa, Bolt bending stress at shank
    F_bolt[i] = sig_tens[i] * (bolt_shank_area/1000**2) # N, Bolt tension load
    M_bolt[i] = sig_bend[i] * ( bolt_Wsh / 1000**3) # N.m, Bolt bending load
    sig_tens_As[i] = F_bolt[i] / (bolt_stress_area/1000**2) # Pa, Axial stress at stress area
    sig_bend_As[i] = M_bolt[i] / (bolt_Ws/1000**3) # Pa, Bending stress at stress area
    delta_sig_total[i] = ((F_bolt[i] - F_bolt[i][0])/ (bolt_stress_area/1000**2)) + ((M_bolt[i] - M_bolt[i][0])/ (bolt_Ws/1000**3))
    delta_sig_tens[i] = ((F_bolt[i] - F_bolt[i][0])/ (bolt_stress_area/1000**2))
    delta_sig_bend[i] = ((M_bolt[i] - M_bolt[i][0])/ (bolt_Ws/1000**3))

	#Schmidt-Neuper Calculation
SN_Fb_diff = SN_Fb - SN_Fb[(SN_Fz == -(abs(SN_Fz).min())) + (SN_Fz == (abs(SN_Fz).min()))]

    # Plot the bolt load at bolt station 1
plt.style.use(['science', 'no-latex'])

times = np.array([0, 0.2, 0.4, 0.6, 0.8, 1.0])
    # Tower load
F_z_bolt = times * 942000 # N, Bolt load for individual bolt, over steps

if imperf_def == "Random":

    fig, ax = plt.subplots(1,1)
    fig.set_size_inches(6, 2.5)
    ax.grid()
    fig.text(0.5, 0.02, 'Tower load per bolt, $F_z$ [kN]', ha = 'center', va='center')
    
    for i in range(0, job_num):
        ax.plot((F_z_bolt)/1000, (probeData["Run_"+str(i)][(len(probe_locs)-1)/2]['Load_steSIEQ_NOEU_Vector:0'].resample('0.2S').interpolate()* (bolt_shank_area/1000**2))/1000, color = "k", alpha=Imperf_height[i], label = str(round(Imperf_height[i], 2)) + " mm")

    ax.plot(SN_Fz/1000, SN_Fb/1000, color = "orange", label = "Schmidt-Neuper")
    ax.set_xlim([0,1000])

    ax.set_ylim([2150,2550])
    ax.set_ylabel('Bolt axial load, $F_s$ [kN]')

    hands, labs = ax.get_legend_handles_labels()
    hands = [x for _, x in sorted(zip(labs, hands), key=lambda pair: pair[0])]
    labs.sort()

    ax.legend(handles = hands, labels = labs, loc='upper left', ncol=3)
    
    fig.savefig(run_filepath+'/Bolt_load.pdf')

else:

    fig, ax = plt.subplots(1,2)
    fig.set_size_inches(6, 2.5)
    ax[0].grid()
    ax[1].grid()
    ax[1].yaxis.set_label_position("right")
    ax[1].yaxis.tick_right()
    fig.text(0.5, 0.02, 'Tower load per bolt, $F_z$ [kN]', ha = 'center', va='center')
    markers = ['o', 's', '^','*', '*', '>', 'D', '<']
    labels = ['$0\degree$ / 0 mm','$30\degree$ / 2 mm', '$120\degree$ / 2 mm', ' ' ,'$25\degree$ / 2 mm', 'Schmidt-Neuper']

    for i in range(0, job_num):
        # Get quasi-static timesteps for load interpolation

            # Tower load
        F_z_bolt = times * 942000 # N, Bolt load for individual bolt, over steps
        
        ax[0].scatter((F_z_bolt)/1000, (probeData["Run_"+str(i)][(len(probe_locs)-1)/2]['Load_steSIEQ_NOEU_Vector:0'].resample('0.2S').interpolate()* (bolt_shank_area/1000**2))/1000, marker=markers[i], s=10, clip_on=False, label = labels[i])
        ax[1].scatter((F_z_bolt)/1000, delta_sig_tens[i]/10**6, marker=markers[i], s=10, clip_on=False, label = labels[i])

    ax[0].plot(SN_Fz/1000, SN_Fb/1000, clip_on= True, label = labels[5], c = 'orange')
    ax[1].plot(SN_Fz/1000, ((SN_Fb_diff)/ (bolt_stress_area/1000**2))/(10**6) , clip_on= True, label = labels[5], c = 'orange')

    ax[0].set_xlim([0,1000])
    ax[1].set_xlim([0,1000])

    ax[0].set_ylim([2100,2800])
    ax[1].set_ylim([0,200])

    ax[0].set_ylabel('Bolt axial load, $F_s$ [kN]')
    ax[1].set_ylabel('Axial stress range, $\Delta\sigma_s$ [MPa]')

    ax[0].legend(loc='upper left')
    ax[1].legend(loc='upper left')

    fig.savefig(run_filepath+'/Bolt_load.pdf')

# endregion
    
# region Fatigue comparison
    # import Fatigue module
from Post_pro import Fatigue

    # Set detail category
det_cat = '36*'

    # Calculate damage
Stress_bins_FEA = {}
Stress_bin_SN = {}

Dam_FEA = {}

        # Schmidt-Neuper
Stress_bin_SN[0] = np.ones((len(F_z_bolt),1))
Stress_bin_SN[1] = np.interp(F_z_bolt, SN_Fz, ((SN_Fb_diff)/ (bolt_stress_area/1000**2))/(10**6))

Dam_cumul, Dam_SN, Dam_blocks_single  = Fatigue.EN1993_fat(Stress_bin_SN, bolt_shank_dia, det_cat)

        # FEA
for i in range(0, job_num):
    Stress_bins_FEA[0] = np.ones((len(delta_sig_tens[i][:]),1))
    Stress_bins_FEA[1] = np.array(delta_sig_tens[i][:]/10**6)

    Dam_cumul, Dam_blocks, Dam_blocks_single  = Fatigue.EN1993_fat(Stress_bins_FEA, bolt_shank_dia, det_cat)
    Dam_FEA[i] = Dam_blocks

    # Plot damages
fig, ax = plt.subplots(1,1)
fig.set_size_inches(6, 2.5)
ax.grid()
ax.set_xlabel("Fatigue damage relative to Schmidt-Neuper method, " r"$\frac{D}{D_{SN}}$" " [-]")
Dam_F = []

for i in [0,1,2,4]:
    Dam_F.append(Dam_FEA[i][-1][0])

#Dam_F = [x for _, x in sorted(zip(Imperf_height, Dam_F), key=lambda pair: pair[0])]

labels = [str(Imperf_angle[i]) + "$\degree$ / "+ str(round(Imperf_height[i],2)) + " mm" for i in [0,1,2,4]]
labels = [x for _, x in sorted(zip(Dam_F, labels), key=lambda pair: pair[0])]
Dam_F.sort()

rects = ax.barh(range(0,4), Dam_F/Dam_SN[-1], color = 'grey')#, alpha=Imperf_height[i])

ax.set_yticks(range(len(labels)))
ax.set_yticklabels(labels)
ax.set_xscale("log")
ax.set_ylabel("Imperfection sizes")
ax.axvline(x=1, color='red', linestyle='--', linewidth=2) # Adjust bar

dam_rand = np.loadtxt('/data/Torque_2024/Runs/Axial_runs_30deg_random/Rand_Damages.csv')

ax.axvline(x=dam_rand.min()/Dam_SN[-1], color='blue', linestyle='-.', linewidth=2) # Adjust bar
ax.axvline(x=dam_rand.max()/Dam_SN[-1], color='blue', linestyle='-.', linewidth=2) # Adjust bar

ax.set_xlim([0, 10])

fig.savefig(run_filepath+'/Fat_dam.pdf')

# endregion