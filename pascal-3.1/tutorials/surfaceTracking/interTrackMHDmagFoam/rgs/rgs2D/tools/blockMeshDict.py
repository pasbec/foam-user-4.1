#!/usr/bin/python
# July 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

from foamTools.blockMeshDict import *
import math as m

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

geo_scale = 1e-3

# Background mesh
geo_x1 = -150.0
geo_x2 =  150.0

geo_dy =    1.25 # 5 mm / (2*2)
geo_y1 = geo_dy/2.0
geo_y2 = geo_y1 + geo_dy

geo_z1 =  -60.0
geo_z2 =  140.0

# --------------------------------------------------------------------------- #
# --- blockMeshDict --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

pHeader(geo_scale)

if pSubDict("vertices"):

    # Front
    pVert(  0, [geo_x1, geo_y1, geo_z1])
    pVert(  1, [geo_x2, geo_y1, geo_z1])
    pVert(  2, [geo_x1, geo_y1, geo_z2])
    pVert(  3, [geo_x2, geo_y1, geo_z2])

    # Back
    pVert(  4, [geo_x1, geo_y2, geo_z1])
    pVert(  5, [geo_x2, geo_y2, geo_z1])
    pVert(  6, [geo_x1, geo_y2, geo_z2])
    pVert(  7, [geo_x2, geo_y2, geo_z2])

if pSubDict("blocks"):

    pBlock( 0, [0, 1, 5, 4, 2, 3, 7, 6], [60, 1, 40], z="background")

pSubDict("edges")

if pSubDict("boundary"):

    if pBoundarySubDict("front", "empty"):

        pFace([0, 1, 3, 2])

    if pBoundarySubDict("back", "empty"):

        pFace([5, 4, 6, 7])

    if pBoundarySubDict("infinity", "patch"):

        pFace([4, 0, 2, 6])
        pFace([1, 5, 7, 3])
        pFace([4, 5, 1, 0])
        pFace([2, 3, 7, 6])

pSubDict("mergePatchPairs")

###

pFooter()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
