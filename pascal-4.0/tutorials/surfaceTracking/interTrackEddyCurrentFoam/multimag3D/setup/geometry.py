#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

csb = os.path.basename(os.path.realpath(sys.argv[0]))
csd = os.path.dirname(os.path.realpath(sys.argv[0]))
csn = os.path.splitext(csb)[0]

sys.path.append(os.environ['FOAM_USER_TOOLS'] + '/' + 'python')
sys.path.append('/usr/lib/freecad/lib')

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
# --- Geometry -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Create new document
App.newDocument(csn)
App.setActiveDocument(csn)
d = App.activeDocument()

# Inner cylinder sketch
s1 = d.addObject('Sketcher::SketchObject', 'Sketch1')

s1.Label = 'inner'
s1.Placement = Placement(Vector(0.0, 0.0, par.geo_z1), Rotation(0.0, 0.0, 0, 1.0))
s1.addGeometry(Circle(Vector(0.0, 0.0, 0.0), Vector(0.0, 0.0, 1.0), par.geo_r1)) # Element 0
s1.addConstraint(Constraint('Coincident', -1, 1, 0, 3)) # Contstraint 0
s1.addConstraint(Constraint('Radius', 0, par.geo_r1)) # Constraint 1

# Outer cylinder sketch
s2 = d.addObject('Sketcher::SketchObject', 'Sketch2')

s2.Label = 'outer'
s2.Placement = Placement(Vector(0.0, 0.0, par.geo_z1), Rotation(0.0, 0.0, 0, 1.0))
s2.addGeometry(Circle(Vector(0.0, 0.0, 0.0), Vector(0.0, 0.0, 1.0), par.geo_r2)) # Element 0
s2.addConstraint(Constraint('Coincident', -1, 1, 0, 3)) # Contstraint 0
s2.addConstraint(Constraint('Radius', 0, par.geo_r2)) # Constraint 1

# --------------------------------------------------------------------------- #

# Region: region_fluid
region_fluid = d.addObject('Part::Extrusion', 'Fluid')
region_fluid.Label = 'region_fluid'
region_fluid.Base = s1
region_fluid.Dir = (0.0, 0.0, par.geo_z3)
region_fluid.Solid = True
region_fluid.TaperAngle = 0.0

# Region: region_above
region_above = d.addObject('Part::Extrusion', 'Above')
region_above.Label = 'region_above'
region_above.Base = s2
region_above.Dir = (0.0, 0.0, par.geo_z4)
region_above.Solid = True
region_above.TaperAngle = 0.0

# Region: region_below
region_below = d.addObject('Part::Extrusion', 'Below')
region_below.Label = 'region_below'
region_below.Base = s2
region_below.Dir = (0.0, 0.0, par.geo_z0)
region_below.Solid = True
region_below.TaperAngle = 0.0

# Region: buffer
region_buffer = d.addObject("Part::Cut", "Buffer")
region_buffer.Label = 'region_buffer'
region_buffer.Base = region_above
region_buffer.Tool = region_fluid

# Recompute to get solids and their faces
d.recompute()

# Patch: fixedMesh
patch_fixedMesh = d.addObject('Part::Feature', 'FixedMesh')
patch_fixedMesh.Label = 'patch_fixedMesh'
patch_fixedMesh.Shape = Shell([region_buffer.Shape.Face1,
                               region_buffer.Shape.Face2,
                               region_buffer.Shape.Face3])

# Patch: sideWall
patch_sideWall = d.addObject('Part::Feature', 'SideWall')
patch_sideWall.Label = 'patch_sideWall'
patch_sideWall.Shape = Shell([region_fluid.Shape.Face1])

# Patch: bottomWall
patch_bottomWall = d.addObject('Part::Feature', 'BottomWall')
patch_bottomWall.Label = 'patch_bottomWall'
patch_bottomWall.Shape = Shell([region_fluid.Shape.Face2])

# Patch: trackedSurface
patch_trackedSurface = d.addObject('Part::Feature', 'TrackedSurface')
patch_trackedSurface.Label = 'patch_trackedSurface'
patch_trackedSurface.Shape = Shell([region_fluid.Shape.Face3])

# --------------------------------------------------------------------------- #

# Recompute
d.recompute()

# --------------------------------------------------------------------------- #
# --- Export ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

exportObj = [region_fluid, region_buffer,
             patch_fixedMesh, patch_sideWall,
             patch_bottomWall, patch_trackedSurface]

for e in exportObj:

    m = Mesh.Mesh(e.Shape.tessellate(0.1))

    m.transform(scale)

    m.write(par.dir_triSurface + '/' + csn + '_' + e.Label + '.stl')

# --------------------------------------------------------------------------- #
# --- Save ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

# Save document
d.saveAs(par.dir_triSurface + '/' + csn + '.fcstd')

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
