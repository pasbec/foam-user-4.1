#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# June 2016
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

from foamTools.blockMeshDict import blockMeshDict

import math as m

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

geo_scale = 1e-2

# Background mesh
geo_r0 =  1.5
geo_r1 =  3.0
geo_r2 = 12.0

geo_z0 = -6.0
geo_z1 =  0.0
geo_z2 =  3.0
geo_z3 =  6.0
geo_z4 = 12.0

a = 2.0**(-0.5)
f = 1.2

n_scale = 2.0

nr0 = int(m.ceil(n_scale*(geo_r0)*7))
nr1 = int(m.ceil(n_scale*(geo_r1-geo_r0)*10))
nr2 = int(m.ceil(n_scale*((geo_r2-geo_r1)+(geo_r1-geo_r0))/3.0*10))

nz0 = int(m.ceil(n_scale*(geo_z1-geo_z0)*5))
nz1 = int(m.ceil(n_scale*(geo_z2-geo_z1)*10))
nz2 = int(m.ceil(n_scale*(geo_z3-geo_z2)*10))
nz3 = int(m.ceil(n_scale*(geo_z4-geo_z3)*5))

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

d = blockMeshDict("blockMeshDict")

d.vertices.set(  0, [        0.0,          0,     geo_z0])
d.vertices.set(  1, [     geo_r0,          0,     geo_z0])
d.vertices.set(  2, [ f*a*geo_r0, f*a*geo_r0,     geo_z0])
d.vertices.set(  3, [        0.0,     geo_r0,     geo_z0])
d.vertices.set(  4, [     geo_r1,          0,     geo_z0])
d.vertices.set(  5, [   a*geo_r1,   a*geo_r1,     geo_z0])
d.vertices.set(  6, [        0.0,     geo_r1,     geo_z0])
d.vertices.set(  7, [     geo_r2,          0,     geo_z0])
d.vertices.set(  8, [   a*geo_r2,   a*geo_r2,     geo_z0])
d.vertices.set(  9, [        0.0,     geo_r2,     geo_z0])

baseVertices = [ i for i in range(len(d.vertices.labels)) ]

d.vertices.copyTranslate( 10, baseVertices, [0.0, 0.0, geo_z1-geo_z0])
d.vertices.copyTranslate( 20, baseVertices, [0.0, 0.0, geo_z2-geo_z0])
d.vertices.copyTranslate( 30, baseVertices, [0.0, 0.0, geo_z3-geo_z0])
d.vertices.copyTranslate( 40, baseVertices, [0.0, 0.0, geo_z4-geo_z0])

# Blocks

d.blocks.set(  0, [  0,  1,  2,  3, 10, 11, 12, 13], zone="core")
d.blocks.set(  1, [  1,  4,  5,  2, 11, 14, 15, 12], zone="core")
d.blocks.set(  2, [  3,  2,  5,  6, 13, 12, 15, 16], zone="core")

d.blocks.set(  3, [  4,  7,  8,  5, 14, 17, 18, 15], zone="annulus")
d.blocks.set(  4, [  6,  5,  8,  9, 16, 15, 18, 19], zone="annulus")

baseBlocks = [ i for i in range(len(d.blocks.labels)) ]

d.blocks.copyShiftVerticeLabels(  10, baseBlocks,  10)
d.blocks.copyShiftVerticeLabels(  20, baseBlocks,  20)
d.blocks.copyShiftVerticeLabels(  30, baseBlocks,  30)

# Distributions

d.blocks.distribution.set(  0, "x", nr0)
d.blocks.distribution.set(  0, "y", nr0)
d.blocks.distribution.set(  1, "x", nr1)
d.blocks.distribution.set(  3, "x", nr2)

d.blocks.distribution.set(  0, "z", nz0)
d.blocks.distribution.set( 10, "z", nz1)
d.blocks.distribution.set( 20, "z", nz2)
d.blocks.distribution.set( 30, "z", nz3)

# Gradings

d.blocks.grading.set( 1, [ 0.5, 1.0, 1.0])
d.blocks.grading.set( 3, [ 7.5, 1.0, 1.0])

d.blocks.grading.set( 0, [1.0, 1.0,  0.13333333333333])
d.blocks.grading.set(10, [1.0, 1.0,  2.0])
d.blocks.grading.set(20, [1.0, 1.0,  0.5])
d.blocks.grading.set(30, [1.0, 1.0,  7.5])

# Boundary faces

d.boundaryFaces.set(  0, "mirror_x",  0, "x-")
d.boundaryFaces.set(  1, "mirror_x",  2, "x-")
d.boundaryFaces.set(  2, "mirror_x",  4, "x-")
d.boundaryFaces.set( 10, "mirror_x", 10, "x-")
d.boundaryFaces.set( 11, "mirror_x", 12, "x-")
d.boundaryFaces.set( 12, "mirror_x", 14, "x-")
d.boundaryFaces.set( 20, "mirror_x", 20, "x-")
d.boundaryFaces.set( 21, "mirror_x", 22, "x-")
d.boundaryFaces.set( 22, "mirror_x", 24, "x-")
d.boundaryFaces.set( 30, "mirror_x", 30, "x-")
d.boundaryFaces.set( 31, "mirror_x", 32, "x-")
d.boundaryFaces.set( 32, "mirror_x", 34, "x-")

d.boundaryFaces.set(  3, "mirror_y",  0, "y-")
d.boundaryFaces.set(  4, "mirror_y",  1, "y-")
d.boundaryFaces.set(  6, "mirror_y",  3, "y-")
d.boundaryFaces.set( 13, "mirror_y", 10, "y-")
d.boundaryFaces.set( 14, "mirror_y", 11, "y-")
d.boundaryFaces.set( 16, "mirror_y", 13, "y-")
d.boundaryFaces.set( 23, "mirror_y", 20, "y-")
d.boundaryFaces.set( 24, "mirror_y", 21, "y-")
d.boundaryFaces.set( 26, "mirror_y", 23, "y-")
d.boundaryFaces.set( 33, "mirror_y", 30, "y-")
d.boundaryFaces.set( 34, "mirror_y", 31, "y-")
d.boundaryFaces.set( 36, "mirror_y", 33, "y-")

d.boundaryFaces.set( 100, "infinity",  0, "z-")
d.boundaryFaces.set( 101, "infinity",  1, "z-")
d.boundaryFaces.set( 102, "infinity",  2, "z-")
d.boundaryFaces.set( 103, "infinity",  3, "z-")
d.boundaryFaces.set( 104, "infinity",  4, "z-")

d.boundaryFaces.set( 130, "infinity", 30, "z+")
d.boundaryFaces.set( 131, "infinity", 31, "z+")
d.boundaryFaces.set( 132, "infinity", 32, "z+")
d.boundaryFaces.set( 133, "infinity", 33, "z+")
d.boundaryFaces.set( 134, "infinity", 34, "z+")

d.boundaryFaces.set( 200, "infinity",  3, "x+")
d.boundaryFaces.set( 201, "infinity",  4, "y+")
d.boundaryFaces.set( 210, "infinity", 13, "x+")
d.boundaryFaces.set( 211, "infinity", 14, "y+")
d.boundaryFaces.set( 220, "infinity", 23, "x+")
d.boundaryFaces.set( 221, "infinity", 24, "y+")
d.boundaryFaces.set( 230, "infinity", 33, "x+")
d.boundaryFaces.set( 231, "infinity", 34, "y+")

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

    d.arc(  4,  5,  6)
    d.arc(  6,  5,  4)

    d.arc( 14, 15, 16)
    d.arc( 16, 15, 14)

    d.arc( 24, 25, 26)
    d.arc( 26, 25, 24)

    d.arc( 34, 35, 36)
    d.arc( 36, 35, 34)

    d.arc( 44, 45, 46)
    d.arc( 46, 45, 44)

    d.arc(  7,  8,  9)
    d.arc(  9,  8,  7)

    d.arc( 17, 18, 19)
    d.arc( 19, 18, 17)

    d.arc( 27, 28, 29)
    d.arc( 29, 28, 27)

    d.arc( 37, 38, 39)
    d.arc( 39, 38, 37)

    d.arc( 47, 48, 49)
    d.arc( 49, 48, 47)

if d.subDict("boundary"):

    if d.boundarySubDict("mirror_x", "patch"):

        d.boundaryFaces.write()

    if d.boundarySubDict("mirror_y", "patch"):

        d.boundaryFaces.write()

    if d.boundarySubDict("infinity", "patch"):

        d.boundaryFaces.write()

    pass

if d.subDict("mergePatchPairs"):

    pass

# --------------------------------------------------------------------------- #

d.footer()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
