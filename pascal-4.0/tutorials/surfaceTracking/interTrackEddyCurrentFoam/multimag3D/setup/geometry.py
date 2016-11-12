#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

__path__ = os.path.realpath(__file__)
__dir__ = os.path.dirname(__path__)

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

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

scale = Matrix()
scale.scale(par.geo_scale, par.geo_scale, par.geo_scale)

# --------------------------------------------------------------------------- #
# --- Document -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

App.newDocument(__name__)
App.setActiveDocument(__name__)
d = App.activeDocument()

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
# --- Export ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

exportObj = [r["fluid"], r["buffer"],
             p["fixedMesh"], p["sideWall"], p["bottomWall"], p["trackedSurface"]]

for e in exportObj:

    mesh = Mesh.Mesh(e.Shape.tessellate(0.1))
    #mesh = MeshPart.meshFromShape(Shape=e.Shape, MaxLength=10)

    mesh.transform(scale)

    mesh.write(par.dir_triSurface + "/" + __name__ + "_" + e.Label + ".stl")

# --------------------------------------------------------------------------- #
# --- Save ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

# Save document
d.saveAs(par.dir_triSurface + "/" + __name__ + ".fcstd")

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
