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

def V0z(i): return Vector(0.0, blockMeshDict.v[i][1], 0.0)
def Vz0(i): return Vector(blockMeshDict.v[i][0], 0.0, 0.0)
def Vrz(i): return Vector(blockMeshDict.v[i][0], blockMeshDict.v[i][1], 0.0)

V = dict()

V['solid'] = [V0z(0), V0z(3), V0z(6), Vrz(5), Vrz(4), Vrz(1), V0z(0)]
V['fluid'] = [V0z(6), V0z(11), Vrz(10), Vrz(9), Vrz(8), Vrz(7), Vrz(4),
              Vrz(5), V0z(6)]
V['vessel'] = [V0z(12), V0z(0), Vrz(1), Vrz(4), Vrz(7), Vrz(8), Vrz(16),
               Vrz(15), Vrz(14), Vrz(13), V0z(12)]
V['heater'] = [V0z(17), V0z(12), Vrz(13), Vrz(14), Vrz(15), Vrz(16), Vrz(25),
               Vrz(24), Vrz(23), Vrz(22), Vrz(21), Vrz(20), Vrz(19), Vrz(18),
               V0z(17)]
V['free'] = [V0z(11), V0z(29), V0z(46), Vrz(45), Vrz(44), Vrz(43), Vrz(42),
             Vrz(41), Vrz(40), Vrz(39), Vrz(38), Vrz(37), Vrz(36), Vrz(35),
             Vrz(34), Vrz(33), Vrz(32), Vrz(31), V0z(30), V0z(17), Vrz(18),
             Vrz(19), Vrz(20), Vrz(21), Vrz(22), Vrz(23), Vrz(24), Vrz(25),
             Vrz(16), Vrz(8), Vrz(9), Vrz(10), V0z(11)]

s = dict()

for k in V.keys():

    name = 'Sketch' + k.capitalize()
    label = 'sketch_' + k

    s[k] = d.addObject('Sketcher::SketchObject', name)
    s[k].Label = label
    s[k].Placement = Placement(Vector(0.0, 0.0, 0.0),
                               Rotation(Vector(1.0, 0.0, 0.0), 90))

    def addPolyLine(sketch, vectors):

        l = len(vectors) - 1

        for i in range(l):

            sketch.addGeometry(Line(vectors[i], vectors[i+1]))

    addPolyLine(s[k], V[k])

# --------------------------------------------------------------------------- #
# --- Regions --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

r = dict()

for k in V.keys():

    name = 'Region' + k.capitalize()
    label = 'region_' + k

    r[k] = d.addObject('Part::Revolution', name)
    r[k].Label = label
    r[k].Source = s[k]
    r[k].Axis = (0.00,0.00,1.00)
    r[k].Base = (0.00,0.00,0.00)
    r[k].Angle = 360.00
    r[k].Solid = True

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

exportObj = []

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
