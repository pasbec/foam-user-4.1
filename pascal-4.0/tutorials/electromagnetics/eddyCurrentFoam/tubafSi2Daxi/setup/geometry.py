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

v['all']  = [30, 35, 40, 46, 30]
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

cv    = dict()

cv[0] = cs*np.array([par.coil_path['r'] - par.coil_bundle['r']/2.0,
                      par.coils_origin[2] - par.coil_bundle['z']/2.0])
cv[1] = cs*np.array([par.coil_path['r'] + par.coil_bundle['r']/2.0,
                      par.coils_origin[2] - par.coil_bundle['z']/2.0])
cv[2] = cs*np.array([par.coil_path['r'] + par.coil_bundle['r']/2.0,
                      par.coils_origin[2] + par.coil_bundle['z']/2.0])
cv[3] = cs*np.array([par.coil_path['r'] - par.coil_bundle['r']/2.0,
                      par.coils_origin[2] + par.coil_bundle['z']/2.0])
cv[4] = cs*np.array([par.coil_path['r'] - par.coil_bundle['r']/2.0,
                      par.coils_origin[2] - par.coil_bundle['z']/2.0])

s['coil'] = d.addObject('Sketcher::SketchObject', 'SketchCoil')
s['coil'].Label = 'sketch_coil'
s['coil'].Placement = Placement(Vector(0.0, 0.0, 0.0),
                                Rotation(Vector(1.0, 0.0, 0.0), 90))

addPolyLine(s['coil'] , cv.keys(), cv)

# --------------------------------------------------------------------------- #
# --- Regions --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

r = dict()

for k in s.keys():

    name = 'Region' + k.capitalize()
    label = 'region_' + k

    back = d.addObject('Part::Revolution', name + 'Back')
    back.Label = label + '_back'
    back.Source = s[k]
    back.Axis = (0.00,0.00,1.00)
    back.Base = (0.00,0.00,0.00)
    back.Angle = par.mesh_angle/2.0
    back.Solid = True

    front = d.addObject('Part::Revolution', name + 'Front')
    front.Label = label + '_front'
    front.Source = s[k]
    front.Axis = (0.00,0.00,1.00)
    front.Base = (0.00,0.00,0.00)
    front.Angle = -par.mesh_angle/2.0
    front.Solid = True

    r[k] = d.addObject("Part::MultiFuse", name)
    r[k].Label = label
    r[k].Shapes = [front, back]

# --------------------------------------------------------------------------- #

r['coils'] = Draft.makeArray(r['coil'], Vector(0.0, 0.0, cs*par.coils_step),
                             Vector(0.0, 0.0, 0.0), par.coils_n, 1,
                             name="RegionCoils")
d.recompute()
r['coils'].Label = 'region_coils'

r['conductor'] = d.addObject("Part::MultiFuse", "RegionConductor")
r['conductor'].Label = 'region_conductor'
r['conductor'].Shapes = [r['solid'], r['fluid'], r['heater']]

r['space'] = d.addObject("Part::MultiFuse", "RegionSpace")
r['space'].Label = 'region_space'
r['space'].Shapes = [r['vessel'], r['free']]

# --------------------------------------------------------------------------- #

r2D = dict()

for k in s.keys():

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

r2D['coils'] = Draft.makeArray(r2D['coil'], Vector(0.0, 0.0, cs*par.coils_step),
                               Vector(0.0, 0.0, 0.0), par.coils_n, 1,
                               name="Coils2D")
d.recompute()
r2D['coils'].Label = 'region_coils_2D'

# --------------------------------------------------------------------------- #

r3D = dict()

for k in s.keys():

    name = 'Region' + k.capitalize() + '3D'
    label = 'region_' + k + '_3D'

    r3D[k] = d.addObject('Part::Revolution', name)
    r3D[k].Label = label
    r3D[k].Source = s[k]
    r3D[k].Axis = (0.00,0.00,1.00)
    r3D[k].Base = (0.00,0.00,0.00)
    r3D[k].Angle = 360.00
    r3D[k].Solid = True

r3D['conductor'] = d.addObject("Part::MultiFuse", "RegionConductor3D")
r3D['conductor'].Label = 'region_conductor_3D'
r3D['conductor'].Shapes = [r3D['solid'], r3D['fluid'], r3D['heater']]

r3D['space'] = d.addObject("Part::MultiFuse", "RegionSpace3D")
r3D['space'].Label = 'region_space_3D'
r3D['space'].Shapes = [r3D['vessel'], r3D['free']]

# --------------------------------------------------------------------------- #

r3D['coils'] = Draft.makeArray(r3D['coil'], Vector(0.0, 0.0, cs*par.coils_step),
                               Vector(0.0, 0.0, 0.0), par.coils_n, 1,
                               name="Coils3D")
d.recompute()
r3D['coils'].Label = 'region_coils_3D'

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Patches --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

p = dict()

p['front'] = d.addObject('Part::Feature', 'PatchFront')
p['front'].Label = 'patch_front'
p['front'].Shape = Shell([r['all'].Shape.Face4])

p['back'] = d.addObject('Part::Feature', 'PatchBack')
p['back'].Label = 'patch_back'
p['back'].Shape = Shell([r['all'].Shape.Face7])

p['infinity'] = d.addObject('Part::Feature', 'PatchInifnity')
p['infinity'].Label = 'patch_infinity'
p['infinity'].Shape = Shell([r['all'].Shape.Face1, r['all'].Shape.Face2,
                             r['all'].Shape.Face3, r['all'].Shape.Face5,
                             r['all'].Shape.Face6, r['all'].Shape.Face8])

#pff['infinity'] = list()
#pfl['infinity'] = [25, 26, 27, 28, 30, 32, 34, 36, 38, 40, 42, 44, 45, 46, 47, 48,
                   #53, 54, 55, 56, 58, 60, 62, 64, 66, 68, 70, 72, 73, 74, 75, 76]

#for i in pfl['infinity']:

    #exec("pff['infinity'].append(r['space'].Shape.Face" + str(i) + ")")

#p['infinity'] = d.addObject('Part::Feature', 'PatchInifnity')
#p['infinity'].Label = 'patch_infinity'
#p['infinity'].Shape = Shell(pff['infinity'])

# --------------------------------------------------------------------------- #

d.recompute()

# --------------------------------------------------------------------------- #
# --- Export ---------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

exportObj = r.values() + r2D.values() + r3D.values() + p.values()

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
