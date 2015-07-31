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
geo_lx1 = 145.0

#geo_ly0 =  10.0
geo_ly0 =   1.0

geo_lz0 =  60.0
geo_lz1 = 145.0

geo_fz0 =  10.0
geo_fz1 =  30.0

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

d = blockMeshDict()

d.vertices.set( 0, [-geo_lx1, -geo_ly0, -geo_lz0])
d.vertices.set( 1, [-geo_lx0, -geo_ly0, -geo_lz0])
d.vertices.set( 2, [ geo_lx0, -geo_ly0, -geo_lz0])
d.vertices.set( 3, [ geo_lx1, -geo_ly0, -geo_lz0])
d.vertices.set( 4, [-geo_lx1,  geo_ly0, -geo_lz0])
d.vertices.set( 5, [-geo_lx0,  geo_ly0, -geo_lz0])
d.vertices.set( 6, [ geo_lx0,  geo_ly0, -geo_lz0])
d.vertices.set( 7, [ geo_lx1,  geo_ly0, -geo_lz0])

d.vertices.set( 8, [-geo_lx1, -geo_ly0,      0.0])
d.vertices.set( 9, [-geo_lx0, -geo_ly0,      0.0])
d.vertices.set(10, [ geo_lx0, -geo_ly0,      0.0])
d.vertices.set(11, [ geo_lx1, -geo_ly0,      0.0])
d.vertices.set(12, [-geo_lx1,  geo_ly0,      0.0])
d.vertices.set(13, [-geo_lx0,  geo_ly0,      0.0])
d.vertices.set(14, [ geo_lx0,  geo_ly0,      0.0])
d.vertices.set(15, [ geo_lx1,  geo_ly0,      0.0])

d.vertices.set(16, [-geo_lx1, -geo_ly0,  geo_fz0])
d.vertices.set(17, [-geo_lx0, -geo_ly0,  geo_fz0])
d.vertices.set(18, [ geo_lx0, -geo_ly0,  geo_fz1])
d.vertices.set(19, [ geo_lx1, -geo_ly0,  geo_fz1])
d.vertices.set(20, [-geo_lx1,  geo_ly0,  geo_fz0])
d.vertices.set(21, [-geo_lx0,  geo_ly0,  geo_fz0])
d.vertices.set(22, [ geo_lx0,  geo_ly0,  geo_fz1])
d.vertices.set(23, [ geo_lx1,  geo_ly0,  geo_fz1])

d.vertices.set(24, [-geo_lx1, -geo_ly0,  geo_lz1])
d.vertices.set(25, [-geo_lx0, -geo_ly0,  geo_lz1])
d.vertices.set(26, [ geo_lx0, -geo_ly0,  geo_lz1])
d.vertices.set(27, [ geo_lx1, -geo_ly0,  geo_lz1])
d.vertices.set(28, [-geo_lx1,  geo_ly0,  geo_lz1])
d.vertices.set(29, [-geo_lx0,  geo_ly0,  geo_lz1])
d.vertices.set(30, [ geo_lx0,  geo_ly0,  geo_lz1])
d.vertices.set(31, [ geo_lx1,  geo_ly0,  geo_lz1])

# Blocks starting with 1 (but this is arbitrary)
d.blocks.set(1, [  0,  1,  5,  4,  8,  9, 13, 12], zone="region_dynamic")
d.blocks.set(2, [  1,  2,  6,  5,  9, 10, 14, 13], zone="region_dynamic")
d.blocks.set(3, [  2,  3,  7,  6, 10, 11, 15, 14], zone="region_dynamic")

d.blocks.set(4, [  8,  9, 13, 12, 16, 17, 21, 20], zone="region_dynamic")
d.blocks.set(5, [  9, 10, 14, 13, 17, 18, 22, 21], zone="region_fluid")
d.blocks.set(6, [ 10, 11, 15, 14, 18, 19, 23, 22], zone="region_dynamic")

d.blocks.set(7, [ 16, 17, 21, 20, 24, 25, 29, 28], grading=[1,1,20], zone="region_dynamic")
d.blocks.set(8, [ 17, 18, 22, 21, 25, 26, 30, 29], grading=[1,1,20], zone="region_dynamic")
d.blocks.set(9, [ 18, 19, 23, 22, 26, 27, 31, 30], grading=[1,1,20], zone="region_dynamic")


#d.blocks.setDivider(4, "x", 65)
#d.blocks.setDivider(6, "x", 65)
#d.blocks.setDivider(2, "z", 30)
#d.blocks.setDivider(8, "z", 60)
#d.blocks.setDivider(5, [35, 1, 10])

d.boundaryFaces.set("front", [1, 2, 3, 4, 5, 6, 7, 8, 9], "y-")
d.boundaryFaces.set("back", [1, 2, 3, 4, 5, 6, 7, 8, 9], "y+")
d.boundaryFaces.set("infinity", [1, 2, 3], "z-")
d.boundaryFaces.set("infinity", [7, 8, 9], "z+")
d.boundaryFaces.set("infinity", [1, 4, 7], "x-")
d.boundaryFaces.set("infinity", [3, 6, 9], "x+")

#quit()

# --------------------------------------------------------------------------- #
# --- blockMeshDict --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

d.header(geo_scale)

# --------------------------------------------------------------------------- #

if d.subDict("vertices"):

    d.vertices.write()

if d.subDict("blocks"):

    d.blocks.write()

if d.subDict("edges"):

    pass

if d.subDict("boundary"):

    if d.boundarySubDict("front", "empty"):

        d.boundaryFaces.write()

    if d.boundarySubDict("back", "empty"):

        d.stdout.write("// This is just to demonstrate the manual write")
        d.boundaryFaces.write()

    if d.boundarySubDict("infinity", "patch"):

        d.boundaryFaces.write()

if d.subDict("mergePatchPairs"):

    pass

# --------------------------------------------------------------------------- #

d.footer()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
