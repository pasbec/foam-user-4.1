#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
# Pascal Beckstein (p.beckstein@hzdr.de)

from __future__ import nested_scopes, generators, division, absolute_import
from __future__ import with_statement, print_function, unicode_literals

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

__name__
__path__ = os.path.realpath(__file__)
__base__ = os.path.basename(__path__)
__dir__ = os.path.dirname(__path__)
__head__ = os.path.splitext(__base__)[0]

sys.path.append(os.environ["FOAM_USER_TOOLS"] + "/" + "python")
sys.path.append("/usr/lib/freecad/lib")

import math as m
import numpy as np

from foamTools.freecad import addPolyLine, faceShell, exportMeshes

import FreeCAD, Sketcher, Draft, Part, PartDesign, Mesh, MeshPart
from FreeCAD import Units, Placement, Matrix, Vector, Rotation
from Part import Line, Circle, Shell
from Sketcher import Constraint

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

import parameters as par
import blockMeshDict

FreeCAD.newDocument(__head__)
FreeCAD.setActiveDocument(__head__)
d = FreeCAD.activeDocument()

# --------------------------------------------------------------------------- #
# --- Sketches -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

v = dict()

v["all"]  = [30, 35, 40, 46]
v["solid"]  = [0, 3, 6, 5, 4, 1]
v["fluid"]  = [6, 11, 10, 9, 8, 7, 4, 5]
v["vessel"] = [12, 0, 1, 4, 7, 8, 16, 15, 14, 13]
v["heater"] = [17, 12, 13, 14, 15, 16, 25, 24, 23, 22, 21, 20, 19, 18]
v["free"]   = [11, 29, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33,
               32, 31, 30, 17, 18, 19, 20, 21, 22, 23, 24, 25, 16, 8, 9, 10]

s = dict()

for k in v.keys():

    name = "Sketch" + k.capitalize()
    label = "sketch_" + k

    s[k] = d.addObject("Sketcher::SketchObject", name)
    s[k].Label = label
    s[k].Placement = Placement(Vector(0.0, 0.0, 0.0),
                               Rotation(Vector(1.0, 0.0, 0.0), 90))

    addPolyLine(s[k], v[k], blockMeshDict.v)

cs = par.coil_scale/par.geo_scale

cv    = dict()

cv[0] = cs*np.array([par.coil_path["r"] - par.coil_bundle["r"]/2.0,
                      par.coils_origin[2] - par.coil_bundle["z"]/2.0])
cv[1] = cs*np.array([par.coil_path["r"] + par.coil_bundle["r"]/2.0,
                      par.coils_origin[2] - par.coil_bundle["z"]/2.0])
cv[2] = cs*np.array([par.coil_path["r"] + par.coil_bundle["r"]/2.0,
                      par.coils_origin[2] + par.coil_bundle["z"]/2.0])
cv[3] = cs*np.array([par.coil_path["r"] - par.coil_bundle["r"]/2.0,
                      par.coils_origin[2] + par.coil_bundle["z"]/2.0])

s["coil"] = d.addObject("Sketcher::SketchObject", "SketchCoil")
s["coil"].Label = "sketch_coil"
s["coil"].Placement = Placement(Vector(0.0, 0.0, 0.0),
                                Rotation(Vector(1.0, 0.0, 0.0), 90))

addPolyLine(s["coil"] , cv.keys(), cv)

# --------------------------------------------------------------------------- #
# --- Bodies ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

bo = dict()

for k in s.keys():

    name = "Body" + k.capitalize()
    label = "body_" + k

    back = d.addObject("Part::Revolution", name + "Back")
    back.Label = label + "_back"
    back.Source = s[k]
    back.Axis = (0.00,0.00,1.00)
    back.Base = (0.00,0.00,0.00)
    back.Angle = par.mesh_angle/2.0
    back.Solid = True

    front = d.addObject("Part::Revolution", name + "Front")
    front.Label = label + "_front"
    front.Source = s[k]
    front.Axis = (0.00,0.00,1.00)
    front.Base = (0.00,0.00,0.00)
    front.Angle = -par.mesh_angle/2.0
    front.Solid = True

    bo[k] = d.addObject("Part::MultiFuse", name)
    bo[k].Label = label
    bo[k].Shapes = [front, back]

bo["conductor"] = d.addObject("Part::MultiFuse", "BodyConductor")
bo["conductor"].Label = "body_conductor"
bo["conductor"].Shapes = [bo["solid"], bo["fluid"], bo["heater"]]

bo["space"] = d.addObject("Part::MultiFuse", "BodySpace")
bo["space"].Label = "body_space"
bo["space"].Shapes = [bo["vessel"], bo["free"]]

bo["coils"] = Draft.makeArray(bo["coil"], Vector(0.0, 0.0, cs*par.coils_step),
                              Vector(0.0, 0.0, 0.0), par.coils_n, 1,
                              name="BodyCoils")
d.recompute()
bo["coils"].Label = "body_coils"

# --------------------------------------------------------------------------- #

bo2D = dict()

for k in s.keys():

    name = "Body" + k.capitalize() + "2D"
    label = "body_" + k + "_2D"

    back = d.addObject("Part::Extrusion", name + "Back")
    back.Label = label + "_back"
    back.Base = s[k]
    back.Dir = (0.0, par.mesh_thickness/2.0, 0.0)
    back.Solid = True
    back.TaperAngle = 0.0

    front = d.addObject("Part::Extrusion", name + "Front")
    front.Label = label + "_front"
    front.Base = s[k]
    front.Dir = (0.0,-par.mesh_thickness/2.0, 0.0)
    front.Solid = True
    front.TaperAngle = 0.0

    bo2D[k] = d.addObject("Part::MultiFuse", name)
    bo2D[k].Label = label
    bo2D[k].Shapes = [front, back]

bo2D["conductor"] = d.addObject("Part::MultiFuse", "BodyConductor2D")
bo2D["conductor"].Label = "body_conductor_2D"
bo2D["conductor"].Shapes = [bo2D["solid"], bo2D["fluid"], bo2D["heater"]]

bo2D["space"] = d.addObject("Part::MultiFuse", "BodySpace2D")
bo2D["space"].Label = "body_space_2D"
bo2D["space"].Shapes = [bo2D["vessel"], bo2D["free"]]

bo2D["coils"] = Draft.makeArray(bo2D["coil"],
                                Vector(0.0, 0.0, cs*par.coils_step),
                                Vector(0.0, 0.0, 0.0), par.coils_n, 1,
                                name="Coils2D")
d.recompute()
bo2D["coils"].Label = "body_coils_2D"

# --------------------------------------------------------------------------- #

bo3D = dict()

for k in s.keys():

    name = "Body" + k.capitalize() + "3D"
    label = "body_" + k + "_3D"

    bo3D[k] = d.addObject("Part::Revolution", name)
    bo3D[k].Label = label
    bo3D[k].Source = s[k]
    bo3D[k].Axis = (0.00,0.00,1.00)
    bo3D[k].Base = (0.00,0.00,0.00)
    bo3D[k].Angle = 360.00
    bo3D[k].Solid = True

bo3D["conductor"] = d.addObject("Part::MultiFuse", "BodyConductor3D")
bo3D["conductor"].Label = "body_conductor_3D"
bo3D["conductor"].Shapes = [bo3D["solid"], bo3D["fluid"], bo3D["heater"]]

bo3D["space"] = d.addObject("Part::MultiFuse", "BodySpace3D")
bo3D["space"].Label = "body_space_3D"
bo3D["space"].Shapes = [bo3D["vessel"], bo3D["free"]]

bo3D["coils"] = Draft.makeArray(bo3D["coil"], Vector(0.0, 0.0, cs*par.coils_step),
                                Vector(0.0, 0.0, 0.0), par.coils_n, 1,
                                name="Coils3D")
d.recompute()
bo3D["coils"].Label = "body_coils_3D"

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Shells ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

shd = dict()

shd["front"] = [bo["all"], [4]]
shd["back"] = [bo["all"], [7]]
shd["infinity"] = [bo["all"], [1, 2, 3, 5, 6, 8]]
shd["topWall"] = [bo["fluid"], [1, 2, 3, 5, 6, 12]]
shd["sideWall"] = [bo["fluid"], [11, 16]]
shd["cornerWall"] = [bo["fluid"], [10, 15]]
shd["bottomWall"] = [bo["fluid"], [8, 9, 13, 14]]

sh = dict()

for k in shd.keys():

    name = "Shell" + k.capitalize()
    label = "shell_" + k

    sh[k] = d.addObject("Part::Feature", name)
    sh[k].Label = label
    sh[k].Shape = faceShell(shd[k])

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Main ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

bodyObj = bo.values()
shellObj = sh.values()

def main():

    exportObj = bodyObj + bo2D.values() + bo3D.values() + shellObj

    #Export body abnd shell meshes
    exportMeshes(exportObj, par.dir_triSurface, __head__, scale=par.geo_scale)

    # Save document
    d.saveAs(par.dir_triSurface + "/" + __head__ + ".fcstd")

# --------------------------------------------------------------------------- #

if __name__ == "__main__": main()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
