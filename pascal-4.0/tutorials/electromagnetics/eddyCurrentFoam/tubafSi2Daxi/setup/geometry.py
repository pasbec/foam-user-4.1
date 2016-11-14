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

from foamTools.freecad import (addPolyLine, makeFuseBody, makeCutBody,
                               makeExtrudeBody, makeRevolveBody,
                               makeOrthoArrayBody, makePolarArrayBody,
                               makeFaceShell, exportMeshes)

import FreeCAD, Sketcher
from FreeCAD import Units, Placement, Matrix, Vector, Rotation
from Part import Line, Circle

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

    back = makeRevolveBody(d, k + "_back", s[k], angle=par.mesh_angle/2.0)

    front = makeRevolveBody(d, k + "_front", s[k], angle=-par.mesh_angle/2.0)

    bo[k] = makeFuseBody(d, k, [front, back])

bo["conductor"] = makeFuseBody(d, "conductor", [bo["solid"],
                                                bo["fluid"],
                                                bo["heater"]])

bo["space"] = makeFuseBody(d, "space", [bo["vessel"], bo["free"]])

bo["coils"] = makeOrthoArrayBody(d, "coils", bo["coil"],
                                 (0.0, 0.0, cs*par.coils_step), par.coils_n)

# --------------------------------------------------------------------------- #

bo2D = dict()

for k in s.keys():

    back = makeExtrudeBody(d, k + "_2D_back", s[k],
                           (0.0, par.mesh_thickness/2.0, 0.0))

    front = makeExtrudeBody(d, k + "_2D_front", s[k],
                            (0.0, -par.mesh_thickness/2.0, 0.0))

    bo2D[k] = makeFuseBody(d, k + "_2D", [front, back])

bo2D["conductor"] = makeFuseBody(d, "conductor_2D", [bo2D["solid"],
                                                     bo2D["fluid"],
                                                     bo2D["heater"]])

bo2D["space"] = makeFuseBody(d, "space_2D", [bo2D["vessel"], bo2D["free"]])

bo2D["coils"] = makeOrthoArrayBody(d, "coils_2D", bo2D["coil"],
                                 (0.0, 0.0, cs*par.coils_step), par.coils_n)

# --------------------------------------------------------------------------- #

bo3D = dict()

for k in s.keys():

    bo3D[k] = makeRevolveBody(d, k + "_3D", s[k])

bo3D["conductor"] = makeFuseBody(d, "conductor_3D", [bo3D["solid"],
                                                     bo3D["fluid"],
                                                     bo3D["heater"]])

bo3D["space"] = makeFuseBody(d, "space_3D", [bo3D["vessel"], bo3D["free"]])

bo2D["coils"] = makeOrthoArrayBody(d, "coils_3D", bo3D["coil"],
                                 (0.0, 0.0, cs*par.coils_step), par.coils_n)

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Shells ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

shd = dict()

shd["front"] = (bo["all"], [4])
shd["back"] = (bo["all"], [7])
shd["infinity"] = (bo["all"], [1, 2, 3, 5, 6, 8])
shd["topWall"] = (bo["fluid"], [1, 2, 3, 5, 6, 12])
shd["sideWall"] = (bo["fluid"], [11, 16])
shd["cornerWall"] = (bo["fluid"], [10, 15])
shd["bottomWall"] = (bo["fluid"], [8, 9, 13, 14])

sh = dict()

for k in shd.keys():

    sh[k] = makeFaceShell(d, k, shd[k])

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
