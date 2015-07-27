#!/usr/bin/python
import math as m

from blockMeshDictFunctions import *

## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ##

geo_scale = 1e-3

# Background mesh
geo_x1 = -150.0
geo_x2 =  150.0

geo_dy =    1.25 # 5 mm / (2*2)
geo_y1 = geo_dy/2.0
geo_y2 = geo_y1 + geo_dy

geo_z1 =  -60.0
geo_z2 =  140.0

## ************************************************************************* ##

pHeader()

pConvertToMeters(geo_scale)

pSubDict("vertices")

# Front
pVert(  0, geo_x1, geo_y1, geo_z1)
pVert(  1, geo_x2, geo_y1, geo_z1)
pVert(  2, geo_x1, geo_y1, geo_z2)
pVert(  3, geo_x2, geo_y1, geo_z2)

# Back
pVert(  4, geo_x1, geo_y2, geo_z1)
pVert(  5, geo_x2, geo_y2, geo_z1)
pVert(  6, geo_x1, geo_y2, geo_z2)
pVert(  7, geo_x2, geo_y2, geo_z2)

pSubDict("blocks")

print "    hex (0 1 5 4 2 3 7 6) background (60 1 40) simpleGrading (1 1 1)"

pSubDict("edges")

###

pSubDict("boundary")

pBoundarySubDict("front", "empty")

print "            (0 1 3 2)"

pBoundarySubDict("back", "empty")

print "            (5 4 6 7)"

pBoundarySubDict("infinity", "patch")

print "            (4 0 2 6)"
print "            (1 5 7 3)"
print "            (4 5 1 0)"
print "            (2 3 7 6)"

pSubDict("mergePatchPairs")

###

pFooter()

## ************************************************************************* ##
