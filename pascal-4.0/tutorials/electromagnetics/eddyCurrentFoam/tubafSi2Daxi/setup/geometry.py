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

import FreeCAD, Sketcher, Draft, Part, PartDesign, Mesh, MeshPart
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

# --------------------------------------------------------------------------- #

def addPolyLine(s, v, v2D):

    l = len(v) - 1

    def V(i): return Vector(v2D[i][0], v2D[i][1], 0.0)

    for i in range(l):

        s.addGeometry(Line(V(v[i]), V(v[i+1])))

# --------------------------------------------------------------------------- #

s = dict()

for k in v.keys():

    name = 'Sketch' + k.capitalize()
    label = 'sketch_' + k

    s[k] = d.addObject('Sketcher::SketchObject', name)
    s[k].Label = label
    s[k].Placement = Placement(Vector(0.0, 0.0, 0.0),
                               Rotation(Vector(1.0, 0.0, 0.0), 90))

    addPolyLine(s[k], v[k], blockMeshDict.v)

# --------------------------------------------------------------------------- #

cs = par.coil_scale/par.geo_scale

c0v    = dict()

c0v[0] = cs*np.array([par.coil_r - par.coil_bundle['x']/2.0,
                      par.coil_z - par.coil_bundle['y']/2.0])
c0v[1] = cs*np.array([par.coil_r + par.coil_bundle['x']/2.0,
                      par.coil_z - par.coil_bundle['y']/2.0])
c0v[2] = cs*np.array([par.coil_r + par.coil_bundle['x']/2.0,
                      par.coil_z + par.coil_bundle['y']/2.0])
c0v[3] = cs*np.array([par.coil_r - par.coil_bundle['x']/2.0,
                      par.coil_z + par.coil_bundle['y']/2.0])
c0v[4] = cs*np.array([par.coil_r - par.coil_bundle['x']/2.0,
                      par.coil_z - par.coil_bundle['y']/2.0])

c0s = d.addObject('Sketcher::SketchObject', 'SketchCoil0')
c0s.Label = 'sketch_coil0'
c0s.Placement = Placement(Vector(0.0, 0.0, 0.0),
                                Rotation(Vector(1.0, 0.0, 0.0), 90))
addPolyLine(c0s, c0v.keys(), c0v)

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

# --------------------------------------------------------------------------- #

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

c0 = d.addObject('Part::Revolution', 'Coil0')
c0.Label = 'coil0'
c0.Source = c0s
c0.Axis = (0.00,0.00,1.00)
c0.Base = (0.00,0.00,0.00)
c0.Angle = 360.00
c0.Solid = True

c = Draft.makeArray(c0, Vector(0.0, 0.0, cs*par.coil_dn),
                    Vector(0.0, 0.0, 0.0), par.coil_n, 1, name="Coil")
d.recompute()
c.Label = 'coil'

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Patches --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

p = dict()

p['infinity'] = d.addObject('Part::Feature', 'PatchInifnity')
p['infinity'].Label = 'patch_infinity'
p['infinity'].Shape = Shell([r['space'].Shape.Face1, r['space'].Shape.Face2,
                             r['space'].Shape.Face3, r['space'].Shape.Face4,
                             r['space'].Shape.Face5, r['space'].Shape.Face6,
                             r['space'].Shape.Face7, r['space'].Shape.Face8,
                             r['space'].Shape.Face9, r['space'].Shape.Face10,
                             r['space'].Shape.Face11, r['space'].Shape.Face12,
                             r['space'].Shape.Face13, r['space'].Shape.Face14,
                             r['space'].Shape.Face15, r['space'].Shape.Face16 ])

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Export ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

exportObj = r.values() + r2D.values() + [c] + p.values()

for e in exportObj:

    mesh = Mesh.Mesh(e.Shape.tessellate(0.1))
    #mesh = MeshPart.meshFromShape(Shape=e.Shape, MaxLength=10)

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
