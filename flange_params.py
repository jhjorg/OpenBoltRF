#!/usr/bin/env python

# Flange_params.py
# Define parameters for flange under study
# Used for geometry creation, meshing, and simulation to ensure consistency

# region Import packages

import math
import numpy as np

# endregion

# region Flange parameters - define flange parameters
    # Flange
flange_in_dia 	= 6355 # mm
flange_out_dia 	= 7000 # mm
flange_width 	= 86.2 # mm
flange_depth 	= abs((flange_in_dia - flange_out_dia)/2) # mm
flange_height 	= 150 # mm
flange_chamfer 	= 0 # mm
flange_fil_radius 	= 20 # mm
flange_lip		= 106 # mm, flange lip height to meet shell
flange_sym_ang  = 180 # degrees, angle of flange sector to model
flange_arran    = 1 # 1 =  inside flange (shell at outer diameter), 0 = outside flange (shell at inner diameter)

    # Shell
shell_width 		= 86.2 # mm
shell_depth 		= 75 # mm
shell_height 		= flange_out_dia/2 # mm
shell_chamfer 		= flange_chamfer # mm
shell_fillet_rad 	= flange_fil_radius # mm

    # washer
washer_height 	= 10 # mm
washer_dia_in   = 74 # mm
washer_dia_out 	= 125 # mm

    # Bolt
bolt_stress_area = 3460 # mm^2
bolt_stress_dia = (bolt_stress_area*4/math.pi)**0.5 # mm
bolt_shank_dia = 72 # mm
bolt_clamp_length = 2 * (flange_height + washer_height) # mm
bolt_length = 400 # mm
bolt_thread_length = 126 # mm
bolt_shank_length = bolt_length - 150 # mm
bolt_Lg = bolt_length - bolt_thread_length # mm
bolt_head_height = 45 # mm
bolt_dia_flat = 110 # mm
bolt_num = 160 # Number of bolts
bolt_angle = 360/bolt_num # angle of segment with one bolt
bolt_cir_dia = 6675 # mm
bolt_a = (bolt_cir_dia - flange_in_dia)/2 # Flange a dimension, inner dia to bolt circle dia
bolt_b = (flange_out_dia - bolt_cir_dia)/2 # Flange b dimensions, outer dia to bolt circle dia
bolt_sym_num = int(bolt_num*(flange_sym_ang/360)) # Number of bolts in flange sector model
bolt_shank_area = math.pi*0.25*bolt_shank_dia**2 # mm2, Bolt shank area
bolt_Ws       =   (np.pi*0.25)*(bolt_stress_dia/2)**3 # mm3, Bolt bending stress constant - stress area
bolt_Wsh      =   (np.pi*0.25)*(bolt_shank_dia/2)**3 # mm3, Bolt bending stress constant - shank area

    # Area of flange segment surface
if flange_arran == 1:
    flange_area     = ((flange_out_dia**2)-((flange_out_dia-2*shell_depth)**2))*math.pi*0.25*(flange_sym_ang/360)
elif flange_arran == 0:
    flange_area     = (((flange_in_dia+2*shell_depth)**2)-(flange_in_dia**2))*math.pi*0.25*(flange_sym_ang/360)

    # Bolt hole
hole_dia = 78 # mm

    # nut
nut_height = 58 # mm
nut_dia_flat = bolt_dia_flat # mm

    # Point dimensions for face referencing
x_hole = -((bolt_cir_dia/2)*math.cos((bolt_angle/2)*math.pi/180) - (flange_in_dia/2)) # x dimension of hole centre
y_hole = (bolt_cir_dia/2)*math.sin((bolt_angle/2)*math.pi/180) # y dimension of hole centre

if flange_arran == 1:
    x_shell = -(((flange_out_dia-shell_depth)/2)*math.cos((bolt_angle/2)*math.pi/180) - (flange_in_dia/2)) # x dimension of shell centre
    y_shell = ((flange_out_dia-shell_depth)/2)*math.sin((bolt_angle/2)*math.pi/180) # y dimension of shell centre
elif flange_arran == 0:
    x_shell = -(((flange_in_dia+shell_depth)/2)*math.cos((bolt_angle/2)*math.pi/180) - (flange_in_dia/2)) # x dimension of shell centre
    y_shell = ((flange_in_dia+shell_depth)/2)*math.sin((bolt_angle/2)*math.pi/180) # y dimension of shell centre    

x_washer = -((((bolt_cir_dia-((washer_dia_in+washer_dia_out))/2)/2)*math.cos((bolt_angle/2)*math.pi/180) - (flange_in_dia/2))) # x dimension of washer
y_washer = ((((bolt_cir_dia-((washer_dia_in+washer_dia_out))/2)/2)*math.sin((bolt_angle/2)*math.pi/180))) # y dimension of washer

x_hole_sec = {}
y_hole_sec = {}
x_washer_sec = {}
y_washer_sec = {}

for i in range(0,int(bolt_num/2)):
    x_hole_sec[i] = -((bolt_cir_dia/2)*math.cos(((i+1)*bolt_angle-bolt_angle/2)*math.pi/180) - (flange_in_dia/2)) # x dimension of hole centre
    y_hole_sec[i] = (bolt_cir_dia/2)*math.sin(((i+1)*bolt_angle-bolt_angle/2)*math.pi/180) # y dimension of hole centre
    x_washer_sec[i] = -((((bolt_cir_dia-((washer_dia_in+washer_dia_out))/2)/2)*math.cos(((i+1)*bolt_angle-bolt_angle/2)*math.pi/180) - (flange_in_dia/2))) # x dimension of washer
    y_washer_sec[i] = ((((bolt_cir_dia-((washer_dia_in+washer_dia_out))/2)/2)*math.sin(((i+1)*bolt_angle-bolt_angle/2)*math.pi/180))) # y dimension of washer

x_flange = -5 # x dimension of flange contact faces
y_flange = 5 # y dimension of flange contact faces

z_load = flange_height + flange_lip + shell_height # z dimensions of load surface
z_top = flange_height # z dimensions of top flange contact
z_bottom = -flange_height # z dimensions of bottom flange contact
z_thread = (z_bottom) - washer_height - (nut_height/2) # z dimension of thread surface

if flange_arran == 1:
    if bolt_num ==1:
        x_sym_out = -(((flange_out_dia-shell_depth)/2)*math.cos((bolt_angle)*math.pi/180) - (flange_in_dia/2)) # x dimension of out sym face
        y_sym_out = ((flange_out_dia-shell_depth)/2)*math.sin((bolt_angle)*math.pi/180) # y dimension of out sym face
    elif bolt_num >=1:
        x_sym_out = -(((flange_out_dia-shell_depth)/2)*math.cos((flange_sym_ang)*math.pi/180) - (flange_in_dia/2)) # x dimension of out sym face
        y_sym_out = ((flange_out_dia-shell_depth)/2)*math.sin((flange_sym_ang)*math.pi/180) # y dimension of out sym face
elif flange_arran == 0:
    if bolt_num == 1:
        x_sym_out = -(((flange_in_dia+shell_depth)/2)*math.cos((bolt_angle)*math.pi/180) - (flange_in_dia/2)) # x dimension of out sym face
        y_sym_out = ((flange_in_dia+shell_depth)/2)*math.sin((bolt_angle)*math.pi/180) # y dimension of out sym face 
    elif bolt_num >= 1:
        x_sym_out = -(((flange_in_dia+shell_depth)/2)*math.cos((flange_sym_ang)*math.pi/180) - (flange_in_dia/2)) # x dimension of out sym face
        y_sym_out = ((flange_in_dia+shell_depth)/2)*math.sin((flange_sym_ang)*math.pi/180) # y dimension of out sym face 
    
    # Probe locations
Probes = {}

z_mid = -(flange_height+flange_lip+shell_height)/1000 # m, z co-ord of bolt middle

bolt_probes = np.array([0, math.ceil(bolt_sym_num/2)-1, bolt_sym_num-1]) # bolts to probe, 0 datum

bolt_rotations = 180 - (bolt_angle/2 + bolt_probes * bolt_angle) # Angles of bolts to probe

probe_locs = [-1, -(2/3), -(1/3), 0, (1/3), (2/3), 1]

x = []
y = []

for i in range(0,len(bolt_rotations)):
    angle = bolt_rotations[i]
    for j in range(0, len(probe_locs)):
        x = ((0.5*bolt_cir_dia+((bolt_shank_dia/2)*probe_locs[j]))*math.cos(angle*math.pi/180))/1000
        y = ((0.5*bolt_cir_dia+((bolt_shank_dia/2)*probe_locs[j]))*math.sin(angle*math.pi/180))/1000
        Probes[j+(i*7)] = [x, y, z_mid]

    # Load and pre-load
        # Pre-load
F_v = 2180000 # N, Pre-load for bolt
load_switch = 1     # Apply moment load (1) or pressure load (0)

        # Load - Moment application or pressure load
if load_switch == 1:

    M_y = 435*10**6*(bolt_sym_num/bolt_num) # MN.m, FLS

else:
    F_z = 942000*bolt_sym_num # N, Tower shell-load
    sig_F_z     = F_z / (flange_area/1000**2) # Pa, Stress on flange segment area from shell load

    # Stress loads
sig_F_v_b   = -1 * F_v / (math.pi * bolt_stress_dia/1000 * (nut_height+washer_height)/1000) # Pa, Stress from pre-lod on thread area

sig_F_v_n   = -1 * sig_F_v_b # Pa, Stress from pre-load on nut thread

sig_F_v     = F_v / (math.pi * (bolt_stress_dia/1000)**2 / 4) # Pa, Stress across bolt stress area from pre-load 

    # Material definitions
        # Shell material - S355
E_X = 200e9 # Pa, Youngs Modulus

sig_y_X = 345e6 # Pa, Yield stress at 0.2 % plastic strain

sig_U_X = 470e6 # Pa, Ultimate Tensile Stress

epsi_y_X = sig_y_X / E_X # Strain at yield

epsi_U_X = 0.15 # Total strain at UTS

E_T_X = (sig_U_X - sig_y_X) / (epsi_U_X - epsi_y_X) # Plastic modulus

nu_X = 0.3 # Poissons Ratio

rho_X = 7870 # kg/m3, Density

        # Flange material - S355
E_F = 200e9

sig_y_F = 275e6

sig_U_F = 450e6

epsi_y_F = sig_y_F / E_F

epsi_U_F = 0.22

E_T_F = (sig_U_F - sig_y_F) / (epsi_U_F - epsi_y_F)

nu_F = 0.3 # Poissons Ratio

rho_F = 7870 # kg/m3, Density


        # Bolt material
E_G = 210e9

sig_y_G = 940e6

sig_U_G = 1040e6

epsi_y_G = sig_y_G / E_G

epsi_U_G = 0.09

E_T_G = (sig_U_G - sig_y_G) / (epsi_U_G - epsi_y_G)

nu_G = 0.3 # Poissons Ratio

rho_G = 7870 # kg/m3, Density

# endregion
