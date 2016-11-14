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

from foamTools.freecad import (addPolyLine, makeFuseBody, makeFaceShell,
                               exportMeshes)

import FreeCAD, Sketcher, Part, Mesh, MeshPart
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

s = dict()

s["inner"] = d.addObject("Sketcher::SketchObject", "SketchInner")
s["inner"].Label = "sketch_inner"
s["inner"].Placement = Placement(Vector(0.0, 0.0, par.geo_z1),
                                 Rotation(0.0, 0.0, 0, 1.0))
s["inner"].addGeometry(Circle(Vector(0.0, 0.0, 0.0),
                              Vector(0.0, 0.0, 1.0), par.geo_r1))

s["outer"] = d.addObject("Sketcher::SketchObject", "SketchOuter")
s["outer"].Label = "sketch_outer"
s["outer"].Placement = Placement(Vector(0.0, 0.0, par.geo_z1),
                                 Rotation(0.0, 0.0, 0, 1.0))
s["outer"].addGeometry(Circle(Vector(0.0, 0.0, 0.0),
                              Vector(0.0, 0.0, 1.0), par.geo_r2))

# --------------------------------------------------------------------------- #
# --- Bodies ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

bo = dict()

bo["fluid"] = d.addObject("Part::Extrusion", "BodyFluid")
bo["fluid"].Label = "body_fluid"
bo["fluid"].Base = s["inner"]
bo["fluid"].Dir = (0.0, 0.0, par.geo_z3)
bo["fluid"].Solid = True
bo["fluid"].TaperAngle = 0.0

bo["above"] = d.addObject("Part::Extrusion", "BodyAbove")
bo["above"].Label = "body_above"
bo["above"].Base = s["outer"]
bo["above"].Dir = (0.0, 0.0, par.geo_z4)
bo["above"].Solid = True
bo["above"].TaperAngle = 0.0

bo["below"] = d.addObject("Part::Extrusion", "BodyBelow")
bo["below"].Label = "body_below"
bo["below"].Base = s["outer"]
bo["below"].Dir = (0.0, 0.0, par.geo_z0)
bo["below"].Solid = True
bo["below"].TaperAngle = 0.0

bo["buffer"] = d.addObject("Part::Cut", "BodyBuffer")
bo["buffer"].Label = "body_buffer"
bo["buffer"].Base = bo["above"]
bo["buffer"].Tool = bo["fluid"]

bo["all"] = makeFuseBody(d, "all", [bo["above"], bo["below"]])
bo["conductor"] = makeFuseBody(d, "conductor", [bo["fluid"]])
bo["space"] = makeFuseBody(d, "space", [bo["buffer"], bo["below"]])

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Shells ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

shd = dict()

shd["infinity"] = (bo["all"], [1, 2, 3, 4])
shd["fixedMesh"] = (bo["buffer"], [1, 2, 3])
shd["sideWall"] = (bo["fluid"], [1])
shd["bottomWall"] = (bo["fluid"], [2])
shd["trackedSurface"] = (bo["fluid"], [3])

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

    exportObj = bodyObj + shellObj

    #Export body abnd shell meshes
    exportMeshes(exportObj, par.dir_triSurface, __head__, scale=par.geo_scale)

    # Save document
    d.saveAs(par.dir_triSurface + "/" + __head__ + ".fcstd")

# --------------------------------------------------------------------------- #

if __name__ == "__main__": main()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
