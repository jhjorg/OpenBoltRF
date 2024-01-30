#!/usr/bin/env python
# -*- coding: utf-8 -*-

# region mesh_sm.py 
## Description 
## OpenBoltRF mesh_sm.py file, creating .med mesh files for all parts
__author__ = "Jack Jorgensen"
__license__ = "GNU GPLv3"
__version__ = "Per Git commit"
__maintainer__ = "Jack Jorgensen"
__email__ = "jack.jorgensen@research.uwa.edu.au"
__status__ = "Development"
# endregion

# region Import packages
    # Salome Meca packages - executed inside Salome Meca container
import salome
import GEOM
from salome.geom import geomBuilder
import math
import SALOMEDS
import SMESH
from salome.smesh import smeshBuilder
import salome_notebook

    # System packages
from pathlib import Path
import sys
import numpy as np

# endregion

# region Set path
path_pre_folder = sys.argv[4]
# endregion

# region Flange Parameters & project definition

sys.path.append(path_pre_folder)

from flange_params import * # Import flange parameters to allow calc set-up

project_Number = sys.argv[2]
job_name = sys.argv[3]

file_Name = "Flange_" + project_Number
run_filepath = path_pre_folder + "Runs/" + job_name + "/"
print(run_filepath)

# endregion

# region Salome initialise
salome.salome_init()
notebook = salome_notebook.NoteBook()
# endregion

# region Geometry component

    # Create new geometry space
geompy = geomBuilder.New()
O = geompy.MakeVertex(0, 0, 0)
OX = geompy.MakeVectorDXDYDZ(1, 0, 0)
OY = geompy.MakeVectorDXDYDZ(0, 1, 0)
OZ = geompy.MakeVectorDXDYDZ(0, 0, 1)

geompy.addToStudy( O, 'O' )
geompy.addToStudy( OX, 'OX' )
geompy.addToStudy( OY, 'OY' )
geompy.addToStudy( OZ, 'OZ' )

Vertex_rotate1 = geompy.MakeVertex((flange_in_dia/2)/1000, 0, 0)
geompy.addToStudy(Vertex_rotate1, 'Vertex_rotate1' )
Vertex_rotate2 = geompy.MakeVertex((flange_in_dia/2)/1000, 0, (flange_height+flange_lip+shell_height)/1000)
geompy.addToStudy(Vertex_rotate2, 'Vertex_rotate2' )
Vector_1       = geompy.MakeVector(Vertex_rotate1, Vertex_rotate2)
geompy.addToStudy(Vector_1, 'Vector_1' )

    # Import step files
        # Import flanges
print(run_filepath + file_Name + "_Bottom_Flange_dual_mat_multi.step")
if flange_sym_ang != 180:
    Bottom_Flange = geompy.ImportSTEP(run_filepath + file_Name + "_Bottom_Flange_dual_mat_multi.step", False, True)
    Top_Flange = geompy.ImportSTEP(run_filepath + file_Name + "_Top_Flange_dual_mat_multi.step", False, True)
elif flange_sym_ang == 180:
    Bottom_Flange = geompy.ImportSTEP(run_filepath + file_Name + "_Bottom_Flange_dual_mat_multi.step", False, True)
    Top_Flange = geompy.ImportSTEP(run_filepath + file_Name + "_Top_Flange_dual_mat_multi.step", False, True)

    # Import nut and bolts
Nut_Washer = {}
Bolt_assem = {}

for i in range(0,bolt_sym_num):    
        # Importing steps
    Nut_Washer[i] = geompy.ImportSTEP(run_filepath + file_Name + "_Nut_washer_multi.step", False, True)
    Bolt_assem[i] = geompy.ImportSTEP(run_filepath + file_Name + "_Bolt_assem_multi.step", False, True)
        # Adding to study
    geompy.addToStudy( Bolt_assem[i], 'Bolt_assem'+str(i) )
    geompy.addToStudy( Nut_Washer[i], 'Nut_Washer'+str(i) )
        # Moving to correct sector
    geompy.Rotate(Bolt_assem[i], Vector_1, -(bolt_angle*(i))*math.pi/180.0)
    geompy.Rotate(Nut_Washer[i], Vector_1, -(bolt_angle*(i))*math.pi/180.0)

    # Partition flange/shell objects
geomObj_27 = geompy.MakeVertex(0, 0, (flange_lip+flange_height)/1000)
geomObj_28 = geompy.MakeVertex(1+flange_out_dia, 1+flange_out_dia, (flange_lip+flange_height)/1000)
geomObj_29 = geompy.MakeVertex(0, 1+flange_out_dia, (flange_lip+flange_height)/1000)
geomObj_30 = geompy.MakePlaneThreePnt(geomObj_27, geomObj_29, geomObj_28, 5)

geomObj_31 = geompy.MakeVertex(0, 0, 0)
geomObj_32 = geompy.MakeVertex(1+flange_out_dia, 1+flange_out_dia, 0)
geomObj_33 = geompy.MakeVertex(0, 1+flange_out_dia, 0)
geomObj_34 = geompy.MakePlaneThreePnt(geomObj_31, geomObj_32, geomObj_33, 5)

Bottom_Flange_1 = geompy.MakePartition([Bottom_Flange], [geomObj_30], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
[Bottom_shell,Bottom_flange] = geompy.ExtractShapes(Bottom_Flange_1, geompy.ShapeType["SOLID"], True)
Top_Flange_1 = geompy.MakePartition([Top_Flange], [geomObj_30], [], [], geompy.ShapeType["SOLID"], 0, [], 0)
[Top_shell,Top_flange] = geompy.ExtractShapes(Top_Flange_1, geompy.ShapeType["SOLID"], True)
geompy.MirrorByPlane(Bottom_Flange_1, geomObj_34)
TF_Faces = geompy.ExtractShapes(Top_Flange_1, geompy.ShapeType["FACE"], True)
BF_Faces = geompy.ExtractShapes(Bottom_Flange_1, geompy.ShapeType["FACE"], True)

    # Create vertex points for applications of spring supports
        # Flange supports
if flange_arran == 1:
    geomObj_51 = geompy.GetSubShape(Top_Flange_1, [13])
elif flange_arran == 0:
    geomObj_51 = geompy.GetSubShape(Top_Flange_1, [11])
Vertex_TF = geompy.MakeVertexOnCurve(geomObj_51, 0.5, True)

if flange_arran == 1:
    geomObj_52 = geompy.GetSubShape(Bottom_Flange_1, [13])
elif flange_arran == 0:
    geomObj_52 = geompy.GetSubShape(Bottom_Flange_1, [11])
Vertex_BF = geompy.MakeVertexOnCurve(geomObj_52, 0.5, True)

        # Nut and bolt supports
geomObj_nut = {}
geomObj_bolt = {}
Vertex_nut1 = {}
Vertex_nut2 = {}
Vertex_bolt1 = {}
Vertex_bolt2 = {}

for i in range(0,bolt_sym_num): 
    geomObj_nut[i] = geompy.GetSubShape(Nut_Washer[i], [22])
    Vertex_nut1[i] = geompy.MakeVertexOnCurve(geomObj_nut[i], 0.5, True)
    Vertex_nut2[i] = geompy.MakeVertexOnCurve(geomObj_nut[i], 1, True)
    geomObj_bolt[i] = geompy.GetSubShape(Bolt_assem[i], [29])
    Vertex_bolt1[i] = geompy.MakeVertexOnCurve(geomObj_bolt[i], 0.5, True)
    Vertex_bolt2[i] = geompy.MakeVertexOnCurve(geomObj_bolt[i], 1, True)

    # Add parts to study
        # Flanges
geompy.addToStudy( Bottom_Flange_1, 'Bottom_Flange_1' )
geompy.addToStudy( Top_Flange_1, 'Top_Flange_1' )

        # Vertex for spring supports
geompy.addToStudy( Vertex_TF, 'Vertex_TF' )
geompy.addToStudy( Vertex_BF, 'Vertex_BF' )

for i in range(0,bolt_sym_num): 
    geompy.addToStudy( Vertex_nut1[i], 'Vertex_nut1_'+str(i) )
    geompy.addToStudy( Vertex_nut2[i], 'Vertex_nut2_'+str(i) )
    geompy.addToStudy( Vertex_bolt1[i], 'Vertex_bolt1_'+str(i) )
    geompy.addToStudy( Vertex_bolt2[i], 'Vertex_bolt2_'+str(i) )

    # Add vertex references to study for face definition
Vertex_Load     = geompy.MakeVertex(x_shell/1000,y_shell/1000,z_load/1000) # Load application face
Vertex_Fix      = geompy.MakeVertex(x_shell/1000,y_shell/1000,-z_load/1000) # Fixed support face

Vertex_f        = geompy.MakeVertex(x_flange/1000,y_flange/1000,0) # Face of flange contacts

Vertex_sym_tf1    = geompy.MakeVertex(-shell_depth/2/1000,0,flange_height/2/1000) # Faces for sym top flange
Vertex_sym_ts1    = geompy.MakeVertex(-shell_depth/2/1000,0,(z_load-flange_height)/1000) # Faces for sym top shell
Vertex_sym_bs1    = geompy.MakeVertex(-shell_depth/2/1000,0,-(z_load-flange_height)/1000) # Faces for sym bottom shell
Vertex_sym_bf1    = geompy.MakeVertex(-shell_depth/2/1000,0,-flange_height/2/1000) # Faces for sym bottom flange

Vertex_sym_tf2    = geompy.MakeVertex(x_sym_out/1000,y_sym_out/1000,flange_height/2/1000) # Faces for sym top / bottom flange - outer
Vertex_sym_ts2    = geompy.MakeVertex(x_sym_out/1000,y_sym_out/1000,(z_load-flange_height)/1000) # Faces for sym top shell - outer
Vertex_sym_bs2    = geompy.MakeVertex(x_sym_out/1000,y_sym_out/1000,-(z_load-flange_height)/1000) # Faces for sym bottom shell - outer
Vertex_sym_bf2    = geompy.MakeVertex(x_sym_out/1000,y_sym_out/1000,-flange_height/2/1000) # Faces for sym top / bottom flange - outer

Vertex_thr_b = {}
Vertex_f_b = {}
Vertex_f_n = {}

for i in range(0, bolt_sym_num):
    Vertex_thr_b[i]    = geompy.MakeVertex(x_hole_sec[i]/1000,y_hole_sec[i]/1000,z_thread/1000) # Thread face bolt / nut
    Vertex_f_b[i]      = geompy.MakeVertex(x_washer_sec[i]/1000,y_washer_sec[i]/1000,flange_height/1000) # Faces for bolt / flange contact
    Vertex_f_n[i]      = geompy.MakeVertex(x_washer_sec[i]/1000,y_washer_sec[i]/1000,-flange_height/1000) # Faces for nut / flange contact
    geompy.addToStudy(Vertex_thr_b[i], 'Vertex_thr_b'+str(i) )
    geompy.addToStudy(Vertex_f_b[i], 'Vertex_f_b'+str(i) )
    geompy.addToStudy(Vertex_f_n[i], 'Vertex_f_n'+str(i) )

geompy.addToStudy(Vertex_Load, 'Vertex_Load' )
geompy.addToStudy(Vertex_Fix, 'Vertex_Fix' )
geompy.addToStudy(Vertex_f, 'Vertex_f' )
geompy.addToStudy(Vertex_sym_tf1 , 'Vertex_sym_tf1 ' )
geompy.addToStudy(Vertex_sym_ts1, 'Vertex_sym_ts1 ' )
geompy.addToStudy(Vertex_sym_bs1, 'Vertex_sym_bs1' )
geompy.addToStudy(Vertex_sym_bf1, 'Vertex_sym_bf1' )
geompy.addToStudy(Vertex_sym_tf2 , 'Vertex_sym_tf2' )
geompy.addToStudy(Vertex_sym_ts2, 'Vertex_sym_ts2' )
geompy.addToStudy(Vertex_sym_bs2, 'Vertex_sym_bs2' )
geompy.addToStudy(Vertex_sym_bf2, 'Vertex_sym_bf2' )

    # Add faces to bolt assembly
BA_Faces = {}
BA_thr = {}
BA_con = {}

for i in range(0,bolt_sym_num):
    BA_Faces[i] = geompy.ExtractShapes(Bolt_assem[i], geompy.ShapeType["FACE"], True)
    ii = 1
    BA_thr[i] = geompy.SubShapeName(geompy.GetFaceNearPoint(Bolt_assem[i], Vertex_thr_b[i]), Bolt_assem[i]) # Get face for bolt threads
    BA_con[i] = geompy.SubShapeName(geompy.GetFaceNearPoint(Bolt_assem[i], Vertex_f_b[i]), Bolt_assem[i]) # Get face for bolt contact
    for aFace in  BA_Faces[i]:       
        name_i = geompy.SubShapeName(aFace, Bolt_assem[i])
        if (name_i != BA_thr[i]) & ((name_i != BA_con[i])):
            name = name_i
        elif (name_i == BA_thr[i]):
            name = 'Bolt_thread'
        elif (name_i == BA_con[i]):
            name = 'Bolt_con'
        Id_Face = geompy.addToStudyInFather(Bolt_assem[i], aFace, name+str(i))
        ii = ii + 1

    # Add faces to nut washer assembly
NW_Faces = {}
NW_thr = {}
NW_con = {}

for i in range(0,bolt_sym_num):
    NW_Faces[i] = geompy.ExtractShapes(Nut_Washer[i], geompy.ShapeType["FACE"], True)
    ii = 1
    NW_thr[i] = geompy.SubShapeName(geompy.GetFaceNearPoint(Nut_Washer[i], Vertex_thr_b[i]), Nut_Washer[i]) # Get face for nut threads
    NW_con[i] = geompy.SubShapeName(geompy.GetFaceNearPoint(Nut_Washer[i], Vertex_f_n[i]), Nut_Washer[i]) # Get face for nut contact
    for aFace in  NW_Faces[i]:
        name_i = geompy.SubShapeName(aFace, Nut_Washer[i])
        if (name_i != NW_thr[i]) & ((name_i != NW_con[i])):
            name = name_i
        elif (name_i == NW_thr[i]):
            name = 'Nut_thread'
        elif (name_i == NW_con[i]):
            name = 'Nut_con'
        Id_Face = geompy.addToStudyInFather(Nut_Washer[i], aFace, name+str(i))
        ii = ii + 1

    # Add shapes & faces to Bottom Flange
geompy.addToStudyInFather( Bottom_Flange_1, Bottom_flange, 'Bottom_flange' )
geompy.addToStudyInFather( Bottom_Flange_1, Bottom_shell, 'Bottom_shell' )

BF_Faces    = geompy.ExtractShapes(Bottom_Flange_1, geompy.ShapeType["FACE"], True)
ii = 1
BF_Fix      = geompy.SubShapeName(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_Fix), Bottom_Flange_1) # Get face for fixed face
BF_n        = geompy.SubShapeName(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_f_n[0]), Bottom_Flange_1) # Get face for nut contact face
BF_f        = geompy.SubShapeName(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_f), Bottom_Flange_1) # Get face for flange contact face
BF_ssym1    = geompy.SubShapeName(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_sym_bs1), Bottom_Flange_1) # Get face for bottom shell sym 1
BF_fsym1    = geompy.SubShapeName(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_sym_bf1), Bottom_Flange_1) # Get face for bottom shell sym 1
BF_ssym2    = geompy.SubShapeName(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_sym_bs2), Bottom_Flange_1) # Get face for bottom shell sym 2
BF_fsym2    = geompy.SubShapeName(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_sym_bf2), Bottom_Flange_1) # Get face for bottom shell sym 2

for aFace in  BF_Faces:
    name_i = geompy.SubShapeName(aFace, Bottom_Flange_1)
    name = name_i
    if (name_i == BF_Fix):
        name = 'Fixed_face'
    elif (name_i == BF_n):
        name = 'B_Flange_nut_con'
    elif (name_i == BF_f):
        name = 'B_Flange_con'
    elif (name_i == BF_ssym1):
        name = 'B_Flange_ssym1'
    elif (name_i == BF_fsym1):
        name = 'B_Flange_fsym1'
    elif (name_i == BF_ssym2):
        name = 'B_Flange_ssym2'
    elif (name_i == BF_fsym2):
        name = 'B_Flange_fsym2'        
    Id_Face = geompy.addToStudyInFather(Bottom_Flange_1, aFace, name)
    ii = ii + 1

    # Add shapes & faces to Top Flange
geompy.addToStudyInFather( Top_Flange_1, Top_flange, 'Top_flange' )
geompy.addToStudyInFather( Top_Flange_1, Top_shell, 'Top_shell' )

TF_Faces    = geompy.ExtractShapes(Top_Flange_1, geompy.ShapeType["FACE"], True)
ii = 1
TF_Load      = geompy.SubShapeName(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_Load), Top_Flange_1) # Get face for Load face
TF_n        = geompy.SubShapeName(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_f_b[0]), Top_Flange_1) # Get face for bolt contact face
TF_f        = geompy.SubShapeName(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_f), Top_Flange_1) # Get face for flange contact face
TF_ssym1    = geompy.SubShapeName(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_sym_ts1), Top_Flange_1) # Get face for Top shell sym 1
TF_fsym1    = geompy.SubShapeName(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_sym_tf1), Top_Flange_1) # Get face for Top shell sym 1
TF_ssym2    = geompy.SubShapeName(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_sym_ts2), Top_Flange_1) # Get face for Top shell sym 2
TF_fsym2    = geompy.SubShapeName(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_sym_tf2), Top_Flange_1) # Get face for Top shell sym 2

for aFace in  TF_Faces:
    name_i = geompy.SubShapeName(aFace, Top_Flange_1)
    name = name_i
    if (name_i == TF_Load):
        name = 'Load_face'
    elif (name_i == TF_n):
        name = 'T_Flange_bolt_con'
    elif (name_i == TF_f):
        name = 'T_Flange_con'
    elif (name_i == TF_ssym1):
        name = 'T_Flange_ssym1'
    elif (name_i == TF_fsym1):
        name = 'T_Flange_fsym1'
    elif (name_i == TF_ssym2):
        name = 'T_Flange_ssym2'
    elif (name_i == TF_fsym2):
        name = 'T_Flange_fsym2'        
    Id_Face = geompy.addToStudyInFather(Top_Flange_1, aFace, name)
    ii = ii + 1

# endregion

# region Mesh components

    # Create mesh builder
smesh = smeshBuilder.New()

    # Nut washers meshing & Extract parameters for meshing and set to required
Nut_washer_mesh = {}
Nut_netgen = {}
NUT_Parameters = {}

for i in range(0,bolt_sym_num):
    Nut_washer_mesh[i] = smesh.Mesh(Nut_Washer[i])
    Nut_netgen[i] = Nut_washer_mesh[i].Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
    NUT_Parameters[i] = Nut_netgen[i].Parameters()
    NUT_Parameters[i].SetMaxSize( 0.025 )
    NUT_Parameters[i].SetMinSize( 5e-09 )
    NUT_Parameters[i].SetSecondOrder( 1 )
    NUT_Parameters[i].SetOptimize( 1 )
    NUT_Parameters[i].SetFineness( 0 )
    NUT_Parameters[i].SetChordalError( -1 )
    NUT_Parameters[i].SetChordalErrorEnabled( 0 )
    NUT_Parameters[i].SetUseSurfaceCurvature( 1 )
    NUT_Parameters[i].SetFuseEdges( 1 )
    NUT_Parameters[i].SetQuadAllowed( 0 )

    # Bolt assembly mesh - same parameters
Bolt_assem_mesh = {}
Bolt_netgen = {}
BOLT_Parameters = {}

for i in range(0,bolt_sym_num):
    Bolt_assem_mesh[i] = smesh.Mesh(Bolt_assem[i])
    Bolt_netgen[i] = Bolt_assem_mesh[i].Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
    BOLT_Parameters[i] = Bolt_netgen[i].Parameters()
    BOLT_Parameters[i].SetMaxSize( 0.025 )
    BOLT_Parameters[i].SetMinSize( 5e-09 )
    BOLT_Parameters[i].SetSecondOrder( 1 )
    BOLT_Parameters[i].SetOptimize( 1 )
    BOLT_Parameters[i].SetFineness( 0 )
    BOLT_Parameters[i].SetChordalError( -1 )
    BOLT_Parameters[i].SetChordalErrorEnabled( 0 )
    BOLT_Parameters[i].SetUseSurfaceCurvature( 1 )
    BOLT_Parameters[i].SetFuseEdges( 1 )
    BOLT_Parameters[i].SetQuadAllowed( 0 )

    # Bottom flange mesh - same parameters
Bottom_Flange_1_1 = smesh.Mesh(Bottom_Flange_1)
BF_Netgen = Bottom_Flange_1_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
BF_Parameters = BF_Netgen.Parameters()
BF_Parameters.SetMaxSize( 0.075 )
BF_Parameters.SetMinSize( 5e-09 )
BF_Parameters.SetSecondOrder( 1 )
BF_Parameters.SetOptimize( 1 )
BF_Parameters.SetFineness( 0 )
BF_Parameters.SetChordalError( -1 )
BF_Parameters.SetChordalErrorEnabled( 0 )
BF_Parameters.SetUseSurfaceCurvature( 1 )
BF_Parameters.SetFuseEdges( 1 )
BF_Parameters.SetQuadAllowed( 0 )

    # Top flange mesh - same parameters
Top_Flange_1_1 = smesh.Mesh(Top_Flange_1)
TF_Netgen = Top_Flange_1_1.Tetrahedron(algo=smeshBuilder.NETGEN_1D2D3D)
TF_Parameters = TF_Netgen.Parameters()
TF_Parameters.SetMaxSize( 0.075 )
TF_Parameters.SetMinSize( 5e-09 )
TF_Parameters.SetSecondOrder( 1 )
TF_Parameters.SetOptimize( 1 )
TF_Parameters.SetFineness( 0 )
TF_Parameters.SetChordalError( -1 )
TF_Parameters.SetChordalErrorEnabled( 0 )
TF_Parameters.SetUseSurfaceCurvature( 1 )
TF_Parameters.SetFuseEdges( 1 )
TF_Parameters.SetQuadAllowed( 0 )

    # Compute all meshes
isDone = Bottom_Flange_1_1.Compute()
isDone = Top_Flange_1_1.Compute()
for i in range(0,bolt_sym_num):
    isDone = Nut_washer_mesh[i].Compute()
    isDone = Bolt_assem_mesh[i].Compute()

    # Create groups for volumes, faces & nodes
        # Nut / Washer
Nut_washer_vol = {}
Nut_thread_face = {}
Nut_con_face = {}
Nut_washer_nod = {}
Nut_thread_nod = {}
Nut_con_nod = {}


for i in range(0,bolt_sym_num):
    Nut_washer_vol[i] = Nut_washer_mesh[i].GroupOnGeom(Nut_Washer[i],'Nut_Washer_vol'+str(i),SMESH.VOLUME)
    Nut_thread_face[i] = Nut_washer_mesh[i].GroupOnGeom(geompy.GetFaceNearPoint(Nut_Washer[i], Vertex_thr_b[i]),'Nut_thread_face'+str(i),SMESH.FACE)
    Nut_con_face[i] = Nut_washer_mesh[i].GroupOnGeom(geompy.GetFaceNearPoint(Nut_Washer[i], Vertex_f_n[i]),'Nut_con_face'+str(i),SMESH.FACE)
    Nut_washer_nod[i] = Nut_washer_mesh[i].GroupOnGeom(Nut_Washer[i],'Nut_Washer_nod'+str(i),SMESH.NODE)
    Nut_thread_nod[i] = Nut_washer_mesh[i].GroupOnGeom(geompy.GetFaceNearPoint(Nut_Washer[i], Vertex_thr_b[i]),'Nut_thread_nod'+str(i),SMESH.NODE)
    Nut_con_nod[i] = Nut_washer_mesh[i].GroupOnGeom(geompy.GetFaceNearPoint(Nut_Washer[i], Vertex_f_n[i]),'Nut_con_nod'+str(i),SMESH.NODE)

        # Bolt Assembly
Bolt_assem_vol = {}
Bolt_thread_face = {}
Bolt_con_face = {}
Bolt_assem_nod = {}
Bolt_thread_nod = {}
Bolt_con_nod = {}

for i in range(0,bolt_sym_num):
    Bolt_assem_vol[i] = Bolt_assem_mesh[i].GroupOnGeom(Bolt_assem[i],'Bolt_assem_vol'+str(i),SMESH.VOLUME)
    Bolt_thread_face[i] = Bolt_assem_mesh[i].GroupOnGeom(geompy.GetFaceNearPoint(Bolt_assem[i], Vertex_thr_b[i]),'Bolt_thread_face'+str(i),SMESH.FACE)
    Bolt_con_face[i] = Bolt_assem_mesh[i].GroupOnGeom(geompy.GetFaceNearPoint(Bolt_assem[i], Vertex_f_b[i]),'Bolt_con_face'+str(i),SMESH.FACE)
    Bolt_assem_nod[i] = Bolt_assem_mesh[i].GroupOnGeom(Bolt_assem[i],'Bolt_assem_nod'+str(i),SMESH.NODE)
    Bolt_thread_nod[i] = Bolt_assem_mesh[i].GroupOnGeom(geompy.GetFaceNearPoint(Bolt_assem[i], Vertex_thr_b[i]),'Bolt_thread_nod'+str(i),SMESH.NODE)
    Bolt_con_nod[i] = Bolt_assem_mesh[i].GroupOnGeom(geompy.GetFaceNearPoint(Bolt_assem[i], Vertex_f_b[i]),'Bolt_con_nod'+str(i),SMESH.NODE)

        # Bottom Flange
Bottom_Flange_1_2 = Bottom_Flange_1_1.GroupOnGeom(Bottom_Flange_1,'Bottom_Flange_1',SMESH.VOLUME)
Bottom_shell_1 = Bottom_Flange_1_1.GroupOnGeom(Bottom_shell,'Bottom_shell',SMESH.VOLUME)
Bottom_flange_1 = Bottom_Flange_1_1.GroupOnGeom(Bottom_flange,'Bottom_flange',SMESH.VOLUME)
B_Flange_nut_con_1 = Bottom_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_f_n[0]),'B_Flange_nut_con',SMESH.FACE)
B_Flange_ssym1_1 = Bottom_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_sym_bs1),'Sym_1',SMESH.FACE)
B_Flange_con_1 = Bottom_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_f),'B_Flange_con',SMESH.FACE)
B_Flange_fsym1_1 = Bottom_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_sym_bf1),'Sym_1',SMESH.FACE)
B_Flange_ssym2_1 = Bottom_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_sym_bs2),'Sym_2',SMESH.FACE)
BF_Fix = Bottom_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_Fix),'BF_Fix',SMESH.FACE)
B_Flange_fsym2_1 = Bottom_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_sym_bf2),'Sym_2',SMESH.FACE)
Bottom_Flange_1_3 = Bottom_Flange_1_1.GroupOnGeom(Bottom_Flange_1,'Bottom_Flange_1',SMESH.NODE)
Bottom_shell_2 = Bottom_Flange_1_1.GroupOnGeom(Bottom_shell,'Bottom_shell',SMESH.NODE)
Bottom_flange_2 = Bottom_Flange_1_1.GroupOnGeom(Bottom_flange,'Bottom_flange',SMESH.NODE)
B_Flange_nut_con_2 = Bottom_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_f_n[0]),'B_Flange_nut_con',SMESH.NODE)
B_Flange_ssym1_2 = Bottom_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_sym_bs1),'Sym_1',SMESH.NODE)
B_Flange_con_2 = Bottom_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_f),'B_Flange_con',SMESH.NODE)
B_Flange_fsym1_2 = Bottom_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_sym_bf1),'Sym_1',SMESH.NODE)
B_Flange_ssym2_2 = Bottom_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_sym_bs2),'Sym_2',SMESH.NODE)
BF_Fix_2 = Bottom_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_Fix),'BF_Fix',SMESH.NODE)
B_Flange_fsym2_2 = Bottom_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Bottom_Flange_1, Vertex_sym_bf2),'Sym_2',SMESH.NODE)

        # Top Flange
Top_Flange_1_2 = Top_Flange_1_1.GroupOnGeom(Top_Flange_1,'Top_Flange_1',SMESH.VOLUME)
Top_shell_1 = Top_Flange_1_1.GroupOnGeom(Top_shell,'Top_shell',SMESH.VOLUME)
Top_flange_1 = Top_Flange_1_1.GroupOnGeom(Top_flange,'Top_flange',SMESH.VOLUME)
T_Flange_bolt_con_1 = Top_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_f_b[0]),'T_Flange_bolt_con',SMESH.FACE)
T_Flange_ssym1_1 = Top_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_sym_ts1),'Sym_1',SMESH.FACE)
T_Flange_con_1 = Top_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_f),'T_Flange_con',SMESH.FACE)
T_Flange_fsym1_1 = Top_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_sym_tf1),'Sym_1',SMESH.FACE)
T_Flange_ssym2_1 = Top_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_sym_ts2),'Sym_2',SMESH.FACE)
TF_Load = Top_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_Load),'TF_Load',SMESH.FACE)
T_Flange_fsym2_1 = Top_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_sym_tf2),'Sym_2',SMESH.FACE)
Top_Flange_1_3 = Top_Flange_1_1.GroupOnGeom(Top_Flange_1,'Top_Flange_1',SMESH.NODE)
Top_shell_2 = Top_Flange_1_1.GroupOnGeom(Top_shell,'Top_shell',SMESH.NODE)
Top_flange_2 = Top_Flange_1_1.GroupOnGeom(Top_flange,'Top_flange',SMESH.NODE)
T_Flange_bolt_con_2 = Top_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_f_b[0]),'T_Flange_bolt_con',SMESH.NODE)
T_Flange_ssym1_2 = Top_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_sym_ts1),'Sym_1',SMESH.NODE)
T_Flange_con_2 = Top_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_f),'T_Flange_con',SMESH.NODE)
T_Flange_fsym1_2 = Top_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_sym_tf1),'Sym_1',SMESH.NODE)
T_Flange_ssym2_2 = Top_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_sym_ts2),'Sym_2',SMESH.NODE)
TF_Load_2 = Top_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_Load),'TF_Load',SMESH.NODE)
T_Flange_fsym2_2 = Top_Flange_1_1.GroupOnGeom(geompy.GetFaceNearPoint(Top_Flange_1, Vertex_sym_tf2),'Sym_2',SMESH.NODE)

    # Create compound mesh
Meshes = []
Meshes.append(Bottom_Flange_1_1.GetMesh())
Meshes.append(Top_Flange_1_1.GetMesh())

for i in range(0,bolt_sym_num):
    Meshes.append(Nut_washer_mesh[i].GetMesh())
    Meshes.append(Bolt_assem_mesh[i].GetMesh())

Compound_Mesh_1 = smesh.Concatenate( Meshes, 1, 0, 1e-05, False )
smesh.SetName(Compound_Mesh_1, 'Compound_Mesh_1')

    # Add points for spring support application
        # Top flange point
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.NODE,SMESH.FT_LyingOnGeom,SMESH.FT_Undefined,'Vertex_TF',SMESH.FT_Undefined,SMESH.FT_LogicalOR,0.01)
aCriteria.append(aCriterion)
Node_0D = smesh.GetFilterFromCriteria(aCriteria)
Group_TF = Compound_Mesh_1.GroupOnFilter( SMESH.NODE, 'Group_TF', Node_0D )
Node_0D_TF = Group_TF.GetNodeIDs()
elem0d = Compound_Mesh_1.Add0DElement( Node_0D_TF[0] )

        # Bottom flange point
aCriteria = []
aCriterion = smesh.GetCriterion(SMESH.NODE,SMESH.FT_LyingOnGeom,SMESH.FT_Undefined,'Vertex_BF',SMESH.FT_Undefined,SMESH.FT_LogicalOR,0.01)
aCriteria.append(aCriterion)
Node_0D = smesh.GetFilterFromCriteria(aCriteria)
Group_BF = Compound_Mesh_1.GroupOnFilter( SMESH.NODE, 'Group_BF', Node_0D )
Node_0D_BF = Group_BF.GetNodeIDs()
elem0d = Compound_Mesh_1.Add0DElement( Node_0D_BF[0] )

        # Bolt & nut
            # Nut in

for i in range(0,bolt_sym_num):
    aCriteria = []
    aCriterion = smesh.GetCriterion(SMESH.NODE,SMESH.FT_LyingOnGeom,SMESH.FT_Undefined,'Vertex_nut1_'+str(i),SMESH.FT_Undefined,SMESH.FT_LogicalOR,0.0075)
    aCriteria.append(aCriterion)
    Node_0D2 = smesh.GetFilterFromCriteria(aCriteria)
    Group_2 = Compound_Mesh_1.GroupOnFilter( SMESH.NODE, 'Group_2', Node_0D2 )
    Node_0D_ID2 = Group_2.GetNodeIDs()
    elem0d2 = Compound_Mesh_1.Add0DElement( Node_0D_ID2[0] )

            # Nut out

for i in range(0,bolt_sym_num):
    aCriteria = []
    aCriterion = smesh.GetCriterion(SMESH.NODE,SMESH.FT_LyingOnGeom,SMESH.FT_Undefined,'Vertex_nut2_'+str(i),SMESH.FT_Undefined,SMESH.FT_LogicalOR,0.005)
    aCriteria.append(aCriterion)

    Node_0D3 = smesh.GetFilterFromCriteria(aCriteria)

    Group_3 = Compound_Mesh_1.GroupOnFilter( SMESH.NODE, 'Group_2', Node_0D3 )
    Node_0D_ID3 = Group_3.GetNodeIDs()

    elem0d3 = Compound_Mesh_1.Add0DElement( Node_0D_ID3[0] )

            # Bolt in

for i in range(0,bolt_sym_num):
    aCriteria = []
    aCriterion = smesh.GetCriterion(SMESH.NODE,SMESH.FT_LyingOnGeom,SMESH.FT_Undefined,'Vertex_bolt1_'+str(i),SMESH.FT_Undefined,SMESH.FT_LogicalOR,0.005)
    aCriteria.append(aCriterion)

    Node_0D4 = smesh.GetFilterFromCriteria(aCriteria)

    Group_4 = Compound_Mesh_1.GroupOnFilter( SMESH.NODE, 'Group_2', Node_0D4 )
    Node_0D_ID4 = Group_4.GetNodeIDs()

    elem0d4 = Compound_Mesh_1.Add0DElement( Node_0D_ID4[0] )

            # Bolt out

for i in range(0,bolt_sym_num):
    aCriteria = []
    aCriterion = smesh.GetCriterion(SMESH.NODE,SMESH.FT_LyingOnGeom,SMESH.FT_Undefined,'Vertex_bolt2_'+str(i),SMESH.FT_Undefined,SMESH.FT_LogicalOR,0.005)
    aCriteria.append(aCriterion)

    Node_0D5 = smesh.GetFilterFromCriteria(aCriteria)

    Group_5 = Compound_Mesh_1.GroupOnFilter( SMESH.NODE, 'Group_2', Node_0D5 )
    Node_0D_ID5 = Group_5.GetNodeIDs()

    elem0d5 = Compound_Mesh_1.Add0DElement( Node_0D_ID5[0] )

a0D = Compound_Mesh_1.CreateEmptyGroup( SMESH.ELEM0D, 'a0D' )
nbAdd = a0D.AddFrom( Compound_Mesh_1.GetMesh() )

    # Create groups of bolts & nuts to apply material in consistent way
Bolt_nut_mat = []
for i in range(0,bolt_sym_num):
    Bolt_nut_mat.append(Compound_Mesh_1.GetGroupByName('Bolt_assem_vol'+str(i), elemType = SMESH.VOLUME)[0])
    Bolt_nut_mat.append(Compound_Mesh_1.GetGroupByName('Nut_Washer_vol'+str(i), elemType = SMESH.VOLUME)[0])

Bolt_nut = Compound_Mesh_1.UnionListOfGroups(Bolt_nut_mat, 'Bolt_nut')

    # Move mesh to align to moment axes
Compound_Mesh_1.TranslateObject( Compound_Mesh_1, [ -flange_in_dia/2/1000, 0, -(flange_height + flange_lip + shell_height)/1000 ], 0 )

    # Add node at (0,0,0) for moment application
Moment_node = Compound_Mesh_1.CreateEmptyGroup( SMESH.NODE, 'Moment_node' )
Moment_0d = Compound_Mesh_1.CreateEmptyGroup( SMESH.ELEM0D, 'Moment_0d' )
nodeID = Compound_Mesh_1.AddNode( 0, 0, 0 )
nbAdd = Moment_node.Add( [nodeID] )
elem0d = Compound_Mesh_1.Add0DElement( nodeID )
nbAdd = Moment_0d.Add( [elem0d] )

    # Save Mesh
try:
  Compound_Mesh_1.ExportMED(run_filepath + '/Compound_Mesh_'+file_Name+'.med',auto_groups=0,minor=40,overwrite=1,meshPart=None,autoDimension=1)
  pass
except:
  print('ExportMED() failed. Invalid file name?')

# endregion

if salome.sg.hasDesktop():
  salome.sg.updateObjBrowser()
