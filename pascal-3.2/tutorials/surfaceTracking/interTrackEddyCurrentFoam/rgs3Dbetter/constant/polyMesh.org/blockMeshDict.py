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

# Basic

geo_scale      =  1.0e-3

geo_rad        = m.pi / 180.0
geo_cos45      = m.cos(45*geo_rad)
geo_sin45      = m.sin(45*geo_rad)

geo_f001       =  0.25
geo_f112       =  1.0 # 0.25

#geo_inflate_radius = 1.08
geo_inflate_radius = 1.01

# Melt

geo_melt_r     =  6.0
geo_melt_x     = 35.0
geo_melt_y     = 75.0

# Frame

geo_frame_d    = 10.0
geo_frame_x    = geo_melt_x + geo_frame_d
geo_frame_y    = geo_melt_y + geo_frame_d

# Geometry x

geo_x0         = geo_melt_x - geo_melt_r
geo_x01        = geo_melt_x - geo_melt_r/2.0
geo_x001       = geo_x0 + (geo_f001 + (1.0-geo_f001)*geo_cos45)*geo_melt_r/2.0
geo_x011       = geo_x0 + geo_cos45*geo_melt_r*geo_inflate_radius
geo_x1         = geo_melt_x
geo_x112       = geo_x0 + (geo_f112 + (1.0-geo_f112)*geo_cos45)*(geo_melt_r+geo_frame_d/2.0-1.0) # -substrate
geo_x12        = geo_melt_x + geo_frame_d/2.0
geo_x2         = geo_melt_x + geo_frame_d
geo_x3         = 150.0

# Geometry y

geo_y0         = geo_melt_y - geo_melt_r
geo_y01        = geo_melt_y - geo_melt_r/2.0
geo_y001       = geo_y0 + (geo_f001 + (1.0-geo_f001)*geo_sin45)*geo_melt_r/2.0
geo_y011       = geo_y0 + geo_sin45*geo_melt_r*geo_inflate_radius
geo_y1         = geo_melt_y
geo_y112       = geo_y0 + (geo_f112 + (1.0-geo_f112)*geo_sin45)*(geo_melt_r+geo_frame_d/2.0-1.0) # -substrate
geo_y12        = geo_melt_y + geo_frame_d/2.0
geo_y2         = geo_melt_y + geo_frame_d
geo_y3         = 200.0

# Geometry z

geo_z0         = -60.0
geo_z1         =  -5.0
geo_z2         =   0.0
geo_z3         =  20.0
geo_z4         =  50.0
geo_z5         =  90.0
geo_z6         = 140.0

# Grading and Distribution

geo_grade_ref = 2.0
geo_dist_ref = geo_frame_d/2.0

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

d = blockMeshDict("blockMeshDict")

# Vertices

d.vertices.set(  0, [  0.0,      0.0,      geo_z0])
d.vertices.set(  1, [  geo_x0,   0.0,      geo_z0])
d.vertices.set(  2, [  geo_x0,   geo_y0,   geo_z0])
d.vertices.set(  3, [  0.0,      geo_y0,   geo_z0])
# WARNING: Vertice 4 and 5 are missing. Remember this!
d.vertices.set(  6, [  geo_x0,  geo_y01,   geo_z0])
d.vertices.set(  7, [  0.0,     geo_y01,   geo_z0])
d.vertices.set(  8, [ geo_x01,   geo_y0,   geo_z0])
d.vertices.set(  9, [geo_x001, geo_y001,   geo_z0])
d.vertices.set( 10, [ geo_x01,   0.0,      geo_z0])
d.vertices.set( 11, [  geo_x0,   geo_y1,   geo_z0])
d.vertices.set( 12, [  0.0,      geo_y1,   geo_z0])
d.vertices.set( 13, [geo_x011, geo_y011,   geo_z0])
d.vertices.set( 14, [  geo_x1,   geo_y0,   geo_z0])
d.vertices.set( 15, [  geo_x1,   0.0,      geo_z0])
d.vertices.set( 16, [ geo_x01-geo_melt_r/4.0,  geo_y12-1.0,   geo_z0]) # x0 -> x01 # # -substrate
d.vertices.set( 17, [  0.0,     geo_y12-1.0,   geo_z0]) # -substrate
d.vertices.set( 18, [geo_x112, geo_y112,   geo_z0])
d.vertices.set( 19, [ geo_x12-1.0,  geo_y01-geo_melt_r/4.0,   geo_z0]) # y0 -> y01 # # -substrate
d.vertices.set( 20, [ geo_x12-1.0,   0.0,      geo_z0]) # -substrate
d.vertices.set( 21, [  geo_x1-geo_melt_r/2.0,   geo_y2,   geo_z0]) # x0 -> x1 #
d.vertices.set( 22, [  0.0,      geo_y2,   geo_z0])
d.vertices.set( 23, [ geo_x12-geo_melt_r/4.0-1.0+geo_melt_r/4.0,   geo_y2,   geo_z0]) # # -substrate+shift
d.vertices.set( 24, [  geo_x2,  geo_y12-geo_melt_r/4.0-1.0+geo_melt_r/4.0,   geo_z0]) # # -substrate+shift
d.vertices.set( 25, [  geo_x2,   geo_y2,   geo_z0])
d.vertices.set( 26, [  geo_x2,   geo_y1-geo_melt_r/2.0,   geo_z0]) # y0 -> y1 #
d.vertices.set( 27, [  geo_x2,   0.0,      geo_z0])
d.vertices.set( 28, [  geo_x1-geo_melt_r/2.0,   geo_y3,   geo_z0]) # x0 -> x1 #
d.vertices.set( 29, [  0.0,      geo_y3,   geo_z0])
d.vertices.set( 30, [ geo_x12-geo_melt_r/4.0-1.0+geo_melt_r/4.0,   geo_y3,   geo_z0]) # x112 -> geo_x12 # # -substrate+shift
d.vertices.set( 31, [  geo_x2,   geo_y3,   geo_z0])
d.vertices.set( 32, [  geo_x3,   geo_y2,   geo_z0])
d.vertices.set( 33, [  geo_x3,   geo_y3,   geo_z0])
d.vertices.set( 34, [  geo_x3,  geo_y12-geo_melt_r/4.0-1.0+geo_melt_r/4.0,   geo_z0]) # y112 -> geo_y12 # # -substrate+shift
d.vertices.set( 35, [  geo_x3,   geo_y1-geo_melt_r/2.0,   geo_z0]) # y0 -> y1 #
d.vertices.set( 36, [  geo_x3,   0.0,      geo_z0])

baseVertices = [ i for i in range(4) ] + [ i+6 for i in range(31) ]

d.vertices.copyTranslate(100, baseVertices, [0.0, 0.0, geo_z1-geo_z0])
d.vertices.copyTranslate(200, baseVertices, [0.0, 0.0, geo_z2-geo_z0])
d.vertices.copyTranslate(300, baseVertices, [0.0, 0.0, geo_z3-geo_z0])
d.vertices.copyTranslate(400, baseVertices, [0.0, 0.0, geo_z4-geo_z0])
d.vertices.copyTranslate(500, baseVertices, [0.0, 0.0, geo_z5-geo_z0])
d.vertices.copyTranslate(600, baseVertices, [0.0, 0.0, geo_z6-geo_z0])

# Blocks

d.blocks.set(  0, [  0,   1,   2,   3, 100, 101, 102, 103], zone="fluid")
d.blocks.set(  1, [  3,   2,   6,   7, 103, 102, 106, 107], zone="fluid")
d.blocks.set(  2, [  2,   8,   9,   6, 102, 108, 109, 106], zone="fluid")
d.blocks.set(  3, [  1,  10,   8,   2, 101, 110, 108, 102], zone="fluid")
d.blocks.set(  4, [  7,   6,  11,  12, 107, 106, 111, 112], zone="fluid")
d.blocks.set(  5, [  6,   9,  13,  11, 106, 109, 113, 111], zone="fluid")
d.blocks.set(  6, [  8,  14,  13,   9, 108, 114, 113, 109], zone="fluid")
d.blocks.set(  7, [ 10,  15,  14,   8, 110, 115, 114, 108], zone="fluid")
d.blocks.set(  8, [ 12,  11,  16,  17, 112, 111, 116, 117], zone="frame")
d.blocks.set(  9, [ 11,  13,  18,  16, 111, 113, 118, 116], zone="frame")
d.blocks.set( 10, [ 14,  19,  18,  13, 114, 119, 118, 113], zone="frame")
d.blocks.set( 11, [ 15,  20,  19,  14, 115, 120, 119, 114], zone="frame")
d.blocks.set( 12, [ 17,  16,  21,  22, 117, 116, 121, 122], zone="frame")
d.blocks.set( 13, [ 16,  18,  23,  21, 116, 118, 123, 121], zone="frame")
d.blocks.set( 14, [ 18,  24,  25,  23, 118, 124, 125, 123], zone="frame")
d.blocks.set( 15, [ 19,  26,  24,  18, 119, 126, 124, 118], zone="frame")
d.blocks.set( 16, [ 20,  27,  26,  19, 120, 127, 126, 119], zone="frame")
d.blocks.set( 17, [ 22,  21,  28,  29, 122, 121, 128, 129], zone="background")
d.blocks.set( 18, [ 21,  23,  30,  28, 121, 123, 130, 128], zone="background")
d.blocks.set( 19, [ 23,  25,  31,  30, 123, 125, 131, 130], zone="background")
d.blocks.set( 20, [ 25,  32,  33,  31, 125, 132, 133, 131], zone="background")
d.blocks.set( 21, [ 24,  34,  32,  25, 124, 134, 132, 125], zone="background")
d.blocks.set( 22, [ 26,  35,  34,  24, 126, 135, 134, 124], zone="background")
d.blocks.set( 23, [ 27,  36,  35,  26, 127, 136, 135, 126], zone="background")

baseBlocks = [ i for i in range(24) ]

d.blocks.copyShiftVerticeLabels(100, baseBlocks, 100)
d.blocks.copyShiftVerticeLabels(200, baseBlocks, 200)
d.blocks.copyShiftVerticeLabels(300, baseBlocks, 300)
d.blocks.copyShiftVerticeLabels(400, baseBlocks, 400)
d.blocks.copyShiftVerticeLabels(500, baseBlocks, 500)

# Gradings

d.blocks.grading.set( 20, [geo_grade_ref, geo_grade_ref, 1.0/geo_grade_ref])
d.blocks.grading.set(520, [geo_grade_ref, geo_grade_ref,     geo_grade_ref])

# Distributions

d.blocks.distribution.set(  0, "x", max(1, int(m.floor(geo_melt_x/geo_dist_ref))))
d.blocks.distribution.set( 20, "x", max(1, int(m.floor((geo_x3-geo_x2)/geo_dist_ref/geo_grade_ref))))

d.blocks.distribution.set(  0, "y", max(1, int(m.floor(geo_melt_y/geo_dist_ref))))
d.blocks.distribution.set( 20, "y", max(1, int(m.floor((geo_y3-geo_y2)/geo_dist_ref/geo_grade_ref))))

d.blocks.distribution.set(  0, "z", max(1, int(m.floor((geo_z1-geo_z0)/geo_dist_ref/geo_grade_ref))))
d.blocks.distribution.set(100, "z", max(1, int(m.floor((geo_z2-geo_z1)/geo_dist_ref))))
d.blocks.distribution.set(200, "z", max(1, int(m.floor((geo_z3-geo_z2)/geo_dist_ref))))
d.blocks.distribution.set(300, "z", max(1, int(m.floor((geo_z4-geo_z3)/geo_dist_ref))))
d.blocks.distribution.set(400, "z", max(1, int(m.floor((geo_z5-geo_z4)/geo_dist_ref))))
d.blocks.distribution.set(500, "z", max(1, int(m.floor((geo_z6-geo_z5)/geo_dist_ref/geo_grade_ref))))

## Boundary faces

d.boundaryFaces.set(0, "mirror_x", [0, 1, 4, 8, 12, 17], "x-")
d.boundaryFaces.set(3, "mirror_y", [0, 3, 7, 11, 16, 23], "y-")

d.boundaryFaces.set(1, "infinity", [20, 21, 22, 23], "x+")
d.boundaryFaces.set(4, "infinity", [17, 18, 19, 20], "y+")
d.boundaryFaces.set(5, "infinity", [ i for i in range(24) ], "z-")
d.boundaryFaces.set(6, "infinity", [ i+500 for i in range(24) ], "z+")

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

    if d.boundarySubDict("mirror_x", "patch"):

        d.boundaryFaces.write()

    if d.boundarySubDict("mirror_y", "patch"):

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
