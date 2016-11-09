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

import math as m
import numpy as np

import FreeCAD, Sketcher, Part, Mesh, MeshPart
from FreeCAD import Units, Placement, Matrix, Vector, Rotation
from Part import Line, Circle, Shell
from Sketcher import Constraint

import blockMeshDict

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

App.newDocument(csn)
App.setActiveDocument(csn)
d = App.activeDocument()

# --------------------------------------------------------------------------- #
# --- Sketches -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

v = dict()

v['solid']  = [0, 3, 6, 5, 4, 1, 0]
v['fluid']  = [6, 11, 10, 9, 8, 7, 4, 5, 6]
v['vessel'] = [12, 0, 1, 4, 7, 8, 16, 15, 14, 13, 12]
v['heater'] = [17, 12, 13, 14, 15, 16, 25, 24, 23, 22, 21, 20, 19, 18, 17]
v['free']   = [11, 29, 46, 45, 44, 43, 42, 41, 40, 39, 38, 37, 36, 35, 34, 33,
               32, 31, 30, 17, 18, 19, 20, 21, 22, 23, 24, 25, 16, 8, 9, 10, 11]

s = dict()

for k in v.keys():

    name = 'Sketch' + k.capitalize()
    label = 'sketch_' + k

    s[k] = d.addObject('Sketcher::SketchObject', name)
    s[k].Label = label
    s[k].Placement = Placement(Vector(0.0, 0.0, 0.0),
                               Rotation(Vector(1.0, 0.0, 0.0), 90))

    def addPolyLine(s, v, bmdv2D):

        l = len(v) - 1

        def V(i): return Vector(bmdv2D[i][0], bmdv2D[i][1], 0.0)

        for i in range(l):

            s.addGeometry(Line(V(v[i]), V(v[i+1])))

    addPolyLine(s[k], v[k], blockMeshDict.v)

# --------------------------------------------------------------------------- #
# --- Regions --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

r = dict()

for k in v.keys():

    name = 'Region' + k.capitalize()
    label = 'region_' + k

    r[k] = d.addObject('Part::Revolution', name)
    r[k].Label = label
    r[k].Source = s[k]
    r[k].Axis = (0.00,0.00,1.00)
    r[k].Base = (0.00,0.00,0.00)
    r[k].Angle = 360.00
    r[k].Solid = True

r['conductor'] = d.addObject("Part::MultiFuse", "RegionConductor")
r['conductor'].Label = 'region_conductor'
r['conductor'].Shapes = [r['solid'], r['fluid'], r['heater']]

r['space'] = d.addObject("Part::MultiFuse", "RegionSpace")
r['space'].Label = 'region_space'
r['space'].Shapes = [r['vessel'], r['free']]

r2D = dict()

for k in v.keys():

    name = 'Region' + k.capitalize() + '2D'
    label = 'region_' + k + '_2D'

    back = d.addObject('Part::Extrusion', name + 'Back')
    back.Label = label + '_back'
    back.Base = s[k]
    back.Dir = (0.0, par.mesh_thickness/2.0, 0.0)
    back.Solid = True
    back.TaperAngle = 0.0

    front = d.addObject('Part::Extrusion', name + 'Front')
    front.Label = label + '_front'
    front.Base = s[k]
    front.Dir = (0.0,-par.mesh_thickness/2.0, 0.0)
    front.Solid = True
    front.TaperAngle = 0.0

    r2D[k] = d.addObject("Part::MultiFuse", name)
    r2D[k].Label = label
    r2D[k].Shapes = [front, back]

r2D['conductor'] = d.addObject("Part::MultiFuse", "RegionConductor2D")
r2D['conductor'].Label = 'region_conductor_2D'
r2D['conductor'].Shapes = [r2D['solid'], r2D['fluid'], r2D['heater']]

r2D['space'] = d.addObject("Part::MultiFuse", "RegionSpace2D")
r2D['space'].Label = 'region_space_2D'
r2D['space'].Shapes = [r2D['vessel'], r2D['free']]

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Patches --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

p = dict()

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Export ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

exportObj = r.values() + r2D.values()

for e in exportObj:

    mesh = Mesh.Mesh(e.Shape.tessellate(0.1))

    mesh.transform(scale)

    mesh.write(par.dir_triSurface + '/' + csn + '_' + e.Label + '.stl')

# --------------------------------------------------------------------------- #
# --- Save ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

# Save document
d.saveAs(par.dir_triSurface + '/' + csn + '.fcstd')

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
