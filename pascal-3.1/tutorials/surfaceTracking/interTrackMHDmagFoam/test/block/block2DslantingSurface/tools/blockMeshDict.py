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

geo_lx0 =  35.0
geo_lx1 = 150.0

geo_ly0 =   0.5

geo_lz0 =  60.0
geo_lz1 = 140.0

geo_fz0 =  10.0
geo_fz1 =  30.0

d = cBlockMeshDict()

# --------------------------------------------------------------------------- #
# --- blockMeshDict --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

d.header(geo_scale)

# --------------------------------------------------------------------------- #

if d.subDict("vertices"):

    d.vertices.set( 9, [-geo_lx0,  geo_ly0,     0.0])
    d.vertices.set(10, [ geo_lx0,  geo_ly0,     0.0])
    d.vertices.set(13, [-geo_lx0, -geo_ly0,     0.0])
    d.vertices.set(14, [ geo_lx0, -geo_ly0,     0.0])
    d.vertices.set(17, [-geo_lx0,  geo_ly0, geo_fz0])
    d.vertices.set(18, [ geo_lx0,  geo_ly0, geo_fz1])
    d.vertices.set(21, [-geo_lx0, -geo_ly0, geo_fz0])
    d.vertices.set(22, [ geo_lx0, -geo_ly0, geo_fz1])

    d.vertices.write()

if d.subDict("blocks"):

    d.blocks.set(0, [ 9, 10, 14, 13, 17, 18, 22, 21], [70, 1, 20], zone="fluid")

    d.blocks.write()

d.subDict("edges")

if d.subDict("boundary"):

    if d.boundarySubDict("front", "empty") and False:

         printFace([0, 1, 3, 2])

    if d.boundarySubDict("back", "empty") and False:

         printFace([5, 4, 6, 7])

    if d.boundarySubDict("infinity", "patch") and False:

         printFace([4, 0, 2, 6])
         printFace([1, 5, 7, 3])
         printFace([4, 5, 1, 0])
         printFace([2, 3, 7, 6])

d.subDict("mergePatchPairs")

# --------------------------------------------------------------------------- #

d.footer()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
