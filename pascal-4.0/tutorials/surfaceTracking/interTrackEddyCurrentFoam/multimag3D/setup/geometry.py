#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
# Pascal Beckstein (p.beckstein@hzdr.de)

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
# --- Functions ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def addPolyLine(s, v, v2D):

    l = len(v) - 1

    def V(i): return Vector(v2D[i][0], v2D[i][1], 0.0)

    for i in range(l):

        s.addGeometry(Line(V(v[i]), V(v[i+1])))

    s.addGeometry(Line(V(v[l]), V(v[0])))


def exportMeshes(obj, dir, prefix, tol=0.1, scale=1.0):

    S = Matrix()
    S.scale(scale, scale, scale)

    for o in obj:

        mesh = Mesh.Mesh(o.Shape.tessellate(tol))
        #mesh = MeshPart.meshFromShape(Shape=o.Shape, MaxLength=10)

        mesh.transform(S)

        name = prefix + "_" + o.Label + ".stl"

        mesh.write(dir + "/" + name)

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

bo["all"] = d.addObject("Part::MultiFuse", "BodyAll")
bo["all"].Label = "body_all"
bo["all"].Shapes = [bo["above"], bo["below"]]

bo["conductor"] = d.addObject("Part::MultiFuse", "BodyConductor")
bo["conductor"].Label = "body_conductor"
bo["conductor"].Shapes = [bo["fluid"], bo["fluid"]]

bo["space"] = d.addObject("Part::MultiFuse", "BodySpace")
bo["space"].Label = "body_space"
bo["space"].Shapes = [bo["buffer"], bo["below"]]

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Shells ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

sh = dict()

sh["infinity"] = d.addObject("Part::Feature", "ShellInfinity")
sh["infinity"].Label = "shell_infinity"
sh["infinity"].Shape = Shell([bo['all'].Shape.Face1, bo['all'].Shape.Face2,
                              bo['all'].Shape.Face3, bo['all'].Shape.Face4])

sh["fixedMesh"] = d.addObject("Part::Feature", "ShellFixedMesh")
sh["fixedMesh"].Label = "shell_fixedMesh"
sh["fixedMesh"].Shape = Shell([bo["buffer"].Shape.Face1,
                               bo["buffer"].Shape.Face2,
                               bo["buffer"].Shape.Face3])

sh["sideWall"] = d.addObject("Part::Feature", "ShellSideWall")
sh["sideWall"].Label = "shell_sideWall"
sh["sideWall"].Shape = Shell([bo["fluid"].Shape.Face1])

sh["bottomWall"] = d.addObject("Part::Feature", "ShellBottomWall")
sh["bottomWall"].Label = "shell_bottomWall"
sh["bottomWall"].Shape = Shell([bo["fluid"].Shape.Face2])

sh["trackedSurface"] = d.addObject("Part::Feature", "ShellTrackedSurface")
sh["trackedSurface"].Label = "shell_trackedSurface"
sh["trackedSurface"].Shape = Shell([bo["fluid"].Shape.Face3])

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
