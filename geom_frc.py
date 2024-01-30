#!/usr/bin/env python
# -*- coding: utf-8 -*-

# region geom_frc.py 
## Description 
## OpenBoltRF geom_frc.py file, Creating ring-flange geometry file
__author__ = "Jack Jorgensen"
__license__ = "GNU GPLv3"
__version__ = "Per Git commit"
__maintainer__ = "Jack Jorgensen"
__email__ = "jack.jorgensen@research.uwa.edu.au"
__status__ = "Development"
# endregion

# region  Import packages
	# FreeCAD packages
import FreeCAD
import Part
import Sketcher
import BOPTools.SplitFeatures
import Import
import Draft

	# System packages
import sys

# endregion

 # region Set path
path_pre_folder = sys.argv[5]
# endregion

# region Flange Parameters

sys.path.append(path_pre_folder)

from flange_params import * # Import flange parameters to allow calc set-up
# endregion

# region Create new project & save

project_Number = sys.argv[3]
job_name = sys.argv[4]

file_Name = "Flange_" + project_Number
run_filepath = path_pre_folder + "Runs/" + job_name + "/"

App.newDocument(file_Name)

App.getDocument(file_Name).saveAs(run_filepath + file_Name + ".FCstd")
# endregion

# region Geometry creation
	# Flange
		# Generate sketch
App.activeDocument().addObject('Sketcher::SketchObject', 'Sketch')
App.activeDocument().Sketch.Placement = App.Placement(App.Vector(0.000000, 0.000000, 0.000000), App.Rotation(-0.707107, 0.000000, 0.000000, -0.707107))
App.activeDocument().Sketch.MapMode = "Deactivated"
geoList = []
geoList.append(Part.LineSegment(App.Vector(0.000000,0.000000,0),App.Vector(-flange_depth,0.000000,0)))
geoList.append(Part.LineSegment(App.Vector(-flange_depth,0.000000,0),App.Vector(-flange_depth,flange_height,0)))
geoList.append(Part.LineSegment(App.Vector(-flange_depth,flange_height,0),App.Vector(0,flange_height,0)))
geoList.append(Part.LineSegment(App.Vector(0,flange_height,0),App.Vector(0.000000,0.000000,0)))
App.activeDocument().getObject('Sketch').addGeometry(geoList,False)
conList = []
conList.append(Sketcher.Constraint('Coincident',0,2,1,1))
conList.append(Sketcher.Constraint('Coincident',1,2,2,1))
conList.append(Sketcher.Constraint('Coincident',2,2,3,1))
conList.append(Sketcher.Constraint('Coincident',3,2,0,1))
conList.append(Sketcher.Constraint('Horizontal',0))
conList.append(Sketcher.Constraint('Horizontal',2))
conList.append(Sketcher.Constraint('Vertical',1))
conList.append(Sketcher.Constraint('Vertical',3))
App.activeDocument().getObject('Sketch').addConstraint(conList)

		# Revolve flange
FreeCAD.ActiveDocument.addObject("Part::Revolution","Revolve")
FreeCAD.ActiveDocument.Revolve.Source = FreeCAD.ActiveDocument.Sketch
FreeCAD.ActiveDocument.Revolve.Axis = (0.000000000000000,0.000000000000000,1.000000000000000)
FreeCAD.ActiveDocument.Revolve.Base = (flange_in_dia/2,0.000000000000000,0.000000000000000)
FreeCAD.ActiveDocument.Revolve.Angle = -flange_sym_ang
FreeCAD.ActiveDocument.Revolve.Solid = True
FreeCAD.ActiveDocument.Revolve.AxisLink = None
FreeCAD.ActiveDocument.Revolve.Symmetric = False

	# Shell
		# Generate sketch
App.activeDocument().addObject('Sketcher::SketchObject', 'Sketch1')
App.activeDocument().Sketch1.Placement = App.Placement(App.Vector(0.000000, 0.000000, 0.000000), App.Rotation(-0.707107, 0.000000, 0.000000, -0.707107))
App.activeDocument().Sketch1.MapMode = "Deactivated"
if flange_arran == 1:
	geoList = []
	geoList.append(Part.LineSegment(App.Vector(-flange_depth+shell_depth,flange_height+flange_lip,0),App.Vector(-flange_depth,flange_height+flange_lip,0)))
	geoList.append(Part.LineSegment(App.Vector(-flange_depth,flange_height+flange_lip,0),App.Vector(-flange_depth,shell_height+flange_lip+flange_height,0)))
	geoList.append(Part.LineSegment(App.Vector(-flange_depth,shell_height+flange_lip+flange_height,0),App.Vector(-flange_depth+shell_depth,shell_height+flange_lip+flange_height,0)))
	geoList.append(Part.LineSegment(App.Vector(-flange_depth+shell_depth,shell_height+flange_lip+flange_height,0),App.Vector(-flange_depth+shell_depth,flange_height+flange_lip,0)))
elif flange_arran == 0:
	geoList = []
	geoList.append(Part.LineSegment(App.Vector(-shell_depth,flange_height+flange_lip,0),App.Vector(0,flange_height+flange_lip,0)))
	geoList.append(Part.LineSegment(App.Vector(0,flange_height+flange_lip,0),App.Vector(0,shell_height+flange_lip+flange_height,0)))
	geoList.append(Part.LineSegment(App.Vector(0,shell_height+flange_lip+flange_height,0),App.Vector(-shell_depth,shell_height+flange_lip+flange_height,0)))
	geoList.append(Part.LineSegment(App.Vector(-shell_depth,shell_height+flange_lip+flange_height,0),App.Vector(-shell_depth,flange_height+flange_lip,0)))

App.activeDocument().getObject('Sketch1').addGeometry(geoList,False)
conList = []
conList.append(Sketcher.Constraint('Coincident',0,2,1,1))
conList.append(Sketcher.Constraint('Coincident',1,2,2,1))
conList.append(Sketcher.Constraint('Coincident',2,2,3,1))
conList.append(Sketcher.Constraint('Coincident',3,2,0,1))
conList.append(Sketcher.Constraint('Horizontal',0))
conList.append(Sketcher.Constraint('Horizontal',2))
conList.append(Sketcher.Constraint('Vertical',1))
conList.append(Sketcher.Constraint('Vertical',3))
App.activeDocument().getObject('Sketch1').addConstraint(conList)
App.ActiveDocument.recompute()

		# Revolve shell
FreeCAD.ActiveDocument.addObject("Part::Revolution","Revolve1")
FreeCAD.ActiveDocument.Revolve1.Source = FreeCAD.ActiveDocument.Sketch1
FreeCAD.ActiveDocument.Revolve1.Axis = (0.000000000000000,0.000000000000000,1.000000000000000)
FreeCAD.ActiveDocument.Revolve1.Base = (flange_in_dia/2,0.000000000000000,0.000000000000000)
FreeCAD.ActiveDocument.Revolve1.Angle = -flange_sym_ang
FreeCAD.ActiveDocument.Revolve1.Solid = True
FreeCAD.ActiveDocument.Revolve1.AxisLink = None
FreeCAD.ActiveDocument.Revolve1.Symmetric = False

	# Flange lip
		# Generate sketch
App.activeDocument().addObject('Sketcher::SketchObject', 'Sketch2')
App.activeDocument().Sketch2.Placement = App.Placement(App.Vector(0.000000, 0.000000, 0.000000), App.Rotation(-0.707107, 0.000000, 0.000000, -0.707107))
App.activeDocument().Sketch2.MapMode = "Deactivated"
if flange_arran == 1:
	geoList = []
	geoList.append(Part.LineSegment(App.Vector(-flange_depth+shell_depth,flange_height,0),App.Vector(-flange_depth,flange_height,0)))
	geoList.append(Part.LineSegment(App.Vector(-flange_depth,flange_height,0),App.Vector(-flange_depth,flange_height+flange_lip,0)))
	geoList.append(Part.LineSegment(App.Vector(-flange_depth,flange_height+flange_lip,0),App.Vector(-flange_depth+shell_depth,flange_height+flange_lip,0)))
	geoList.append(Part.LineSegment(App.Vector(-flange_depth+shell_depth,flange_height+flange_lip,0),App.Vector(-flange_depth+shell_depth,flange_height,0)))
if flange_arran == 0:
	geoList = []
	geoList.append(Part.LineSegment(App.Vector(-shell_depth,flange_height,0),App.Vector(0,flange_height,0)))
	geoList.append(Part.LineSegment(App.Vector(0,flange_height,0),App.Vector(0,flange_height+flange_lip,0)))
	geoList.append(Part.LineSegment(App.Vector(0,flange_height+flange_lip,0),App.Vector(-shell_depth,flange_height+flange_lip,0)))
	geoList.append(Part.LineSegment(App.Vector(-shell_depth,flange_height+flange_lip,0),App.Vector(-shell_depth,flange_height,0)))

App.activeDocument().getObject('Sketch2').addGeometry(geoList,False)
conList = []
conList.append(Sketcher.Constraint('Coincident',0,2,1,1))
conList.append(Sketcher.Constraint('Coincident',1,2,2,1))
conList.append(Sketcher.Constraint('Coincident',2,2,3,1))
conList.append(Sketcher.Constraint('Coincident',3,2,0,1))
conList.append(Sketcher.Constraint('Horizontal',0))
conList.append(Sketcher.Constraint('Horizontal',2))
conList.append(Sketcher.Constraint('Vertical',1))
conList.append(Sketcher.Constraint('Vertical',3))
App.activeDocument().getObject('Sketch2').addConstraint(conList)
App.ActiveDocument.recompute()

		# Flange lip revolve
FreeCAD.ActiveDocument.addObject("Part::Revolution","Revolve2")
FreeCAD.ActiveDocument.Revolve2.Source = FreeCAD.ActiveDocument.Sketch2
FreeCAD.ActiveDocument.Revolve2.Axis = (0.000000000000000,0.000000000000000,1.000000000000000)
FreeCAD.ActiveDocument.Revolve2.Base = (flange_in_dia/2,0.000000000000000,0.000000000000000)
FreeCAD.ActiveDocument.Revolve2.Angle = -flange_sym_ang
FreeCAD.ActiveDocument.Revolve2.Solid = True
FreeCAD.ActiveDocument.Revolve2.AxisLink = None
FreeCAD.ActiveDocument.Revolve2.Symmetric = False

App.ActiveDocument.recompute()

	# Fuse flange + shell
App.activeDocument().addObject("Part::MultiFuse","Fusion")
App.activeDocument().Fusion.Shapes = [App.activeDocument().Revolve,App.activeDocument().Revolve1,App.activeDocument().Revolve2]
App.activeDocument().getObject('Fusion').Refine = True
App.ActiveDocument.recompute()

	# Flange weld fillet or chamfer
if flange_arran == 1:
	if flange_chamfer != 0:
		FreeCAD.ActiveDocument.addObject("Part::Chamfer","Chamfer")
		FreeCAD.ActiveDocument.Chamfer.Base = FreeCAD.ActiveDocument.Fusion
		__fillets__ = []
		__fillets__.append((10,shell_chamfer,shell_chamfer))
		FreeCAD.ActiveDocument.Chamfer.Edges = __fillets__
		del __fillets__

		FreeCAD.ActiveDocument.addObject("Part::Fillet","Fillet")
		FreeCAD.ActiveDocument.Fillet.Base = FreeCAD.ActiveDocument.Chamfer
		__fillets__ = []
		__fillets__.append((2,shell_fillet_rad,shell_fillet_rad))
		__fillets__.append((12,shell_fillet_rad,shell_fillet_rad))
		FreeCAD.ActiveDocument.Fillet.Edges = __fillets__
		del __fillets__

	elif flange_chamfer == 0:

		FreeCAD.ActiveDocument.addObject("Part::Fillet","Fillet")
		FreeCAD.ActiveDocument.Fillet.Base = FreeCAD.ActiveDocument.Fusion
		__fillets__ = []
		__fillets__.append((10,shell_fillet_rad,shell_fillet_rad))
		FreeCAD.ActiveDocument.Fillet.Edges = __fillets__
		del __fillets__

elif flange_arran == 0:
	if flange_chamfer != 0:
		FreeCAD.ActiveDocument.addObject("Part::Chamfer","Chamfer")
		FreeCAD.ActiveDocument.Chamfer.Base = FreeCAD.ActiveDocument.Fusion
		__fillets__ = []
		__fillets__.append((9,shell_chamfer,shell_chamfer))
		FreeCAD.ActiveDocument.Chamfer.Edges = __fillets__
		del __fillets__

		FreeCAD.ActiveDocument.addObject("Part::Fillet","Fillet")
		FreeCAD.ActiveDocument.Fillet.Base = FreeCAD.ActiveDocument.Chamfer
		__fillets__ = []
		__fillets__.append((4,shell_fillet_rad,shell_fillet_rad))
		__fillets__.append((12,shell_fillet_rad,shell_fillet_rad))
		FreeCAD.ActiveDocument.Fillet.Edges = __fillets__
		del __fillets__

	elif flange_chamfer == 0:

		FreeCAD.ActiveDocument.addObject("Part::Fillet","Fillet")
		FreeCAD.ActiveDocument.Fillet.Base = FreeCAD.ActiveDocument.Fusion
		__fillets__ = []
		__fillets__.append((9,shell_fillet_rad,shell_fillet_rad))
		FreeCAD.ActiveDocument.Fillet.Edges = __fillets__
		del __fillets__

	# Save document
App.getDocument(file_Name).save()

    # Flange group split - for materials
__shape = Part.getShape(App.activeDocument().getObject('Revolve1'),'',needSubElement=False,refine=False)
App.ActiveDocument.addObject('Part::Feature','Revolve003').Shape=__shape
App.ActiveDocument.recompute()
f = BOPTools.SplitFeatures.makeSlice(name='Slice')
f.Base = [App.ActiveDocument.Fillet, App.ActiveDocument.Revolve003][0]
f.Tools = [App.ActiveDocument.Fillet, App.ActiveDocument.Revolve003][1:]
f.Mode = 'Split'
App.ActiveDocument.recompute()
f.purgeTouched()
App.ActiveDocument.recompute()

App.getDocument(file_Name).save()

	# Bolt Hole
		# Make cut-out
App.ActiveDocument.addObject("Part::Cylinder","bolt_hole")
App.activeDocument().bolt_hole.Radius=hole_dia/2
App.activeDocument().bolt_hole.Height=flange_height
App.activeDocument().bolt_hole.Angle=360.00
App.activeDocument().bolt_hole.Placement=App.Placement(App.Vector(x_hole,y_hole,0.00),App.Rotation(App.Vector(0.00,0.00,1.00),0.00))
App.ActiveDocument.recompute()

        # Array bolt hole
Draft.make_polar_array(App.activeDocument().bolt_hole, bolt_sym_num, -(flange_sym_ang - bolt_angle), center = App.Vector(flange_in_dia/2, 0,0))
#Draft.autogroup(Array)
App.ActiveDocument.recompute()

	# Cut bolt hole out of flange
App.activeDocument().addObject("Part::Cut","Top_Flange")
App.activeDocument().Top_Flange.Base = App.activeDocument().Slice
App.activeDocument().Top_Flange.Tool = App.activeDocument().Array
App.ActiveDocument.recompute()

	# Bottom flange - create and align
__shape = Part.getShape(App.activeDocument().getObject('Top_Flange'),'',needSubElement=False,refine=False)
App.ActiveDocument.addObject('Part::Feature','Bottom_Flange').Shape=__shape
App.activeDocument().Bottom_Flange.Placement=App.Placement(App.Vector(0,0,0), App.Rotation(App.Vector(1,0,0),180), App.Vector(0,0,0)).multiply(App.activeDocument().Bottom_Flange.Placement)
App.activeDocument().Bottom_Flange.Placement=App.Placement(App.Vector(0,0,0), App.Rotation(App.Vector(0,0,1),-flange_sym_ang), App.Vector(flange_in_dia/2,0,0)).multiply(App.activeDocument().Bottom_Flange.Placement)
App.ActiveDocument.recompute()

	# Bolt assembley
    	# Bolt shank
App.ActiveDocument.addObject("Part::Cylinder","bolt_shank")
App.ActiveDocument.bolt_shank.Radius = bolt_shank_dia/2
App.ActiveDocument.bolt_shank.Height = bolt_shank_length
App.ActiveDocument.recompute()

    	# Bolt Head
App.ActiveDocument.addObject("Part::Cylinder","bolt_head")
App.ActiveDocument.bolt_head.Radius = bolt_dia_flat/2
App.ActiveDocument.bolt_head.Height = bolt_head_height
App.ActiveDocument.bolt_head.Placement = App.Placement(App.Vector(0,0,bolt_shank_length),App.Rotation(App.Vector(0,0,1),0))
App.ActiveDocument.recompute()

		# Bolt thread transition cone
App.ActiveDocument.addObject("Part::Cone","Cone")
App.ActiveDocument.Cone.Radius1=bolt_stress_dia/2
App.ActiveDocument.Cone.Radius2=bolt_shank_dia/2
App.ActiveDocument.Cone.Height=bolt_Lg-bolt_shank_length
App.ActiveDocument.Cone.Angle=360.00
App.ActiveDocument.Cone.Placement=App.Placement(App.Vector(0.00,0.00,-(bolt_Lg-bolt_shank_length)),App.Rotation(App.Vector(0.00,0.00,1.00),0.00))

    	# Bolt thread
App.ActiveDocument.addObject("Part::Cylinder","bolt_thread")
App.ActiveDocument.bolt_thread.Radius = bolt_stress_dia/2
App.ActiveDocument.bolt_thread.Height = bolt_length-bolt_Lg
App.ActiveDocument.recompute()
App.ActiveDocument.bolt_thread.Placement=App.Placement(App.Vector(0.00,0.00,-((bolt_Lg-bolt_shank_length)+(bolt_length-bolt_Lg))),App.Rotation(App.Vector(0.00,0.00,1.00),0.00))

    	# Bolt (shank + transition + thread + head)
App.activeDocument().addObject("Part::MultiFuse","Bolt")
App.activeDocument().Bolt.Shapes = [App.activeDocument().bolt_shank,App.activeDocument().bolt_head,App.activeDocument().Cone,App.activeDocument().bolt_thread,]
App.ActiveDocument.recompute()

    # Nut
App.ActiveDocument.addObject("Part::Cylinder","Nut1")
App.ActiveDocument.Nut1.Radius = nut_dia_flat/2
App.ActiveDocument.Nut1.Height = nut_height
App.ActiveDocument.Nut1.Placement = App.Placement(App.Vector(0,0,bolt_shank_length-washer_height*2-nut_height-bolt_clamp_length),App.Rotation(App.Vector(0,0,1),0))
App.ActiveDocument.recompute()

    # Washers
App.ActiveDocument.addObject("Part::Cylinder","Washer1")
App.ActiveDocument.Washer1.Radius = washer_dia_out/2
App.ActiveDocument.Washer1.Height = washer_height
App.ActiveDocument.recompute()
App.ActiveDocument.Washer1.Placement = App.Placement(App.Vector(0,0,bolt_shank_length-washer_height*2-bolt_clamp_length),App.Rotation(App.Vector(0,0,1),0))
__shape2 = Part.getShape(App.activeDocument().getObject('Washer1'),'',needSubElement=False,refine=False)
App.ActiveDocument.addObject('Part::Feature','Washer2').Shape=__shape2
App.ActiveDocument.Washer2.Placement = App.Placement(App.Vector(0,0,bolt_shank_length-washer_height),App.Rotation(App.Vector(0,0,1),0))
App.ActiveDocument.recompute()

	# Position bolt, nut and washers
App.ActiveDocument.Bolt.Placement=App.Placement(App.Vector(x_hole,y_hole,-(bolt_shank_length-flange_height-washer_height)), App.Rotation(App.Vector(0,0,1),0), App.Vector(0,0,0)).multiply(App.ActiveDocument.Bolt.Placement)
App.ActiveDocument.Nut1.Placement=App.Placement(App.Vector(x_hole,y_hole,-flange_height-(bolt_shank_length-bolt_clamp_length-washer_height)), App.Rotation(App.Vector(0,0,1),0), App.Vector(0,0,0)).multiply(App.ActiveDocument.Nut1.Placement)
App.ActiveDocument.Washer1.Placement=App.Placement(App.Vector(x_hole,y_hole,-flange_height-(bolt_shank_length-bolt_clamp_length-washer_height)), App.Rotation(App.Vector(0,0,1),0), App.Vector(0,0,0)).multiply(App.ActiveDocument.Washer1.Placement)
App.ActiveDocument.Washer2.Placement=App.Placement(App.Vector(x_hole,y_hole,-(bolt_shank_length-flange_height-washer_height)), App.Rotation(App.Vector(0,0,1),0), App.Vector(0,0,0)).multiply(App.ActiveDocument.Washer2.Placement)
App.ActiveDocument.recompute()

	# Bolt+washer fuse
App.activeDocument().addObject("Part::Fuse","Bolt_washer")
App.activeDocument().Bolt_washer.Base = App.activeDocument().Bolt
App.activeDocument().Bolt_washer.Tool = App.activeDocument().Washer2

	# Nut + washer fuse & hole cut-out
App.activeDocument().addObject("Part::Fuse","Nut_washer")
App.activeDocument().Nut_washer.Base = App.activeDocument().Nut1
App.activeDocument().Nut_washer.Tool = App.activeDocument().Washer1
App.ActiveDocument.addObject("Part::Cylinder","Nut2")
App.ActiveDocument.Nut2.Radius = bolt_stress_dia/2
App.ActiveDocument.Nut2.Height = nut_height+washer_height
App.ActiveDocument.Nut2.Placement = App.Placement(App.Vector(0,0,bolt_shank_length-washer_height*2-nut_height-bolt_clamp_length),App.Rotation(App.Vector(0,0,1),0))
App.ActiveDocument.recompute()
App.ActiveDocument.Nut2.Placement=App.Placement(App.Vector(x_hole,y_hole,-flange_height-(bolt_shank_length-bolt_clamp_length-washer_height)), App.Rotation(App.Vector(0,0,1),0), App.Vector(0,0,0)).multiply(App.ActiveDocument.Nut2.Placement)
App.ActiveDocument.recompute()
App.activeDocument().addObject("Part::Cut","Nut")
App.activeDocument().Nut.Base = App.activeDocument().Nut_washer
App.activeDocument().Nut.Tool = App.activeDocument().Nut2

	# Split bolt+washer for thread area
f = BOPTools.SplitFeatures.makeSlice(name='Bolt_assem')
f.Base = [App.ActiveDocument.Bolt_washer, App.ActiveDocument.Nut][0]
f.Tools = [App.ActiveDocument.Bolt_washer, App.ActiveDocument.Nut][1:]
f.Mode = 'Split'
App.ActiveDocument.recompute()

App.getDocument(file_Name).save()
# endregion

	# Export Top Flange
	### Begin command Std_Export
__objs__=[]
__objs__.append(FreeCAD.activeDocument().getObject("Top_Flange"))
Import.export(__objs__,run_filepath +"/" + file_Name + "_Top_Flange_dual_mat_multi.step")

del __objs__
	### End command Std_Export

	# Export bottom Flange
	### Begin command Std_Export
__objs__=[]
__objs__.append(FreeCAD.activeDocument().getObject("Bottom_Flange"))
Import.export(__objs__,run_filepath +"/" + file_Name + "_Bottom_Flange_dual_mat_multi.step")

del __objs__
	### End command Std_Export

	# Export Bolt Assembly
	### Begin command Std_Export
__objs__=[]
__objs__.append(FreeCAD.activeDocument().getObject("Bolt_assem"))
Import.export(__objs__,run_filepath +"/" + file_Name + "_Bolt_assem_multi.step")

del __objs__
	### End command Std_Export

	# Export Nut assembly
	### Begin command Std_Export
__objs__=[]
__objs__.append(FreeCAD.activeDocument().getObject("Nut"))
Import.export(__objs__,run_filepath +"/" + file_Name + "_Nut_washer_multi.step")

del __objs__
	### End command Std_Export

