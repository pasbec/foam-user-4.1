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
# --- Regions --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

r = dict()

r["fluid"] = d.addObject("Part::Extrusion", "RegionFluid")
r["fluid"].Label = "region_fluid"
r["fluid"].Base = s["inner"]
r["fluid"].Dir = (0.0, 0.0, par.geo_z3)
r["fluid"].Solid = True
r["fluid"].TaperAngle = 0.0

r["above"] = d.addObject("Part::Extrusion", "RegionAbove")
r["above"].Label = "region_above"
r["above"].Base = s["outer"]
r["above"].Dir = (0.0, 0.0, par.geo_z4)
r["above"].Solid = True
r["above"].TaperAngle = 0.0

r["below"] = d.addObject("Part::Extrusion", "RegionBelow")
r["below"].Label = "region_below"
r["below"].Base = s["outer"]
r["below"].Dir = (0.0, 0.0, par.geo_z0)
r["below"].Solid = True
r["below"].TaperAngle = 0.0

r["buffer"] = d.addObject("Part::Cut", "RegionBuffer")
r["buffer"].Label = "region_buffer"
r["buffer"].Base = r["above"]
r["buffer"].Tool = r["fluid"]

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Patches --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

p = dict()

p["fixedMesh"] = d.addObject("Part::Feature", "PatchFixedMesh")
p["fixedMesh"].Label = "patch_fixedMesh"
p["fixedMesh"].Shape = Shell([r["buffer"].Shape.Face1,
                              r["buffer"].Shape.Face2,
                              r["buffer"].Shape.Face3])

p["sideWall"] = d.addObject("Part::Feature", "PatchSideWall")
p["sideWall"].Label = "patch_sideWall"
p["sideWall"].Shape = Shell([r["fluid"].Shape.Face1])

p["bottomWall"] = d.addObject("Part::Feature", "PatchBottomWall")
p["bottomWall"].Label = "patch_bottomWall"
p["bottomWall"].Shape = Shell([r["fluid"].Shape.Face2])

p["trackedSurface"] = d.addObject("Part::Feature", "PatchTrackedSurface")
p["trackedSurface"].Label = "patch_trackedSurface"
p["trackedSurface"].Shape = Shell([r["fluid"].Shape.Face3])

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Main ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

bodyObj = [r["fluid"], r["buffer"]]
shellObj = [p["fixedMesh"], p["sideWall"], p["bottomWall"], p["trackedSurface"]]

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
