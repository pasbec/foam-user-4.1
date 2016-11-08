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

import math as m
import numpy as np

from foamTools.expansion import expansion_de_e, expansion_n_ds
from foamTools.blockMeshDict import blockMeshDict

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

import parameters as p

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

d = blockMeshDict(p.dir_polyMesh + '/' + 'blockMeshDict')

# Vertices

def vr0(r): return np.array([r, 0.0])
def v0z(z): return np.array([p.mesh_R['axis'], z])
def vrz(r, z): return np.array([r, z])
def vsa(s, a): return s * np.array([m.cos(a), m.sin(a)])
def vct(*args):
    print args
    ct = np.zeros(2)
    for a in args: ct = np.add(ct, a)
    return ct/len(args)

v      = dict()

v[  0] = v0z(0.0)
v[  6] = v0z(p.geo_Z['solid'])
v[  3] = vct(v[0], v[6])
v[  4] = vrz(p.geo_R['solid'], v[6][1])
v[  1] = vct(v[0], v[4])
v[  2] = vct(v[0], v[4], v[6])
v[  5] = vct(v[4], v[6])

v[  7] = vrz(p.geo_R[0], p.geo_Z[5])
v[  8] = vrz(v[7][0], p.geo_Z[6])
v[  9] = vrz(v[4][0], v[8][1])
v[ 10] = vrz(v[5][0], v[8][1])
v[ 11] = v0z(v[8][1])

v[ 12] = v0z(p.geo_Z[1])
v[ 13] = v[1] + p.geo_G
v[ 14] = v[4] + p.geo_G
v[ 15] = vrz(p.geo_R[1], p.geo_Z[4])
v[ 16] = vrz(v[15][0], v[8][1])

v[ 17] = v0z(p.geo_Z[0])
v[ 18] = vrz(v[13][0], v[17][1])
v[ 19] = vrz(v[14][0], v[17][1])
v[ 20] = vrz(v[15][0], v[17][1])
v[ 21] = vrz(p.geo_R[2], v[17][1])
v[ 22] = vrz(v[21][0], v[15][1])
v[ 23] = vrz(v[21][0], v[8][1])
v[ 24] = vrz(v[21][0], p.geo_Z[7])
v[ 25] = vrz(v[15][0], v[24][1])
v[ 26] = vrz(v[7][0], v[24][1])
v[ 27] = vrz(v[9][0], v[24][1])
v[ 28] = vrz(v[5][0], v[24][1])
v[ 29] = v0z(v[28][1])

v[ 30] = v0z(p.mesh_Z['C'] - p.mesh_Z['inf'])
v[ 31] = vrz(v[18][0], v[30][1])
v[ 32] = vrz(v[19][0], v[30][1])
v[ 33] = vrz(v[20][0], v[30][1])
v[ 34] = vrz(v[21][0], v[30][1])
v[ 35] = vrz(p.mesh_R['inf'], v[30][1])
v[ 36] = vrz(v[35][0], v[21][1])
v[ 37] = vrz(v[36][0], v[22][1])
v[ 38] = vrz(v[36][0], v[23][1])
v[ 39] = vrz(v[36][0], v[24][1])
v[ 40] = vrz(v[36][0], p.mesh_Z['C'] + p.mesh_Z['inf'])
v[ 41] = vrz(v[24][0], v[40][1])
v[ 42] = vrz(v[25][0], v[40][1])
v[ 43] = vrz(v[26][0], v[40][1])
v[ 44] = vrz(v[27][0], v[40][1])
v[ 45] = vrz(v[28][0], v[40][1])
v[ 46] = vrz(v[29][0], v[40][1])

n = len(v)

def vset(n, i, v2D):

    k = n + 1

    def x(r): return r * m.cos(p.mesh_phi/2.0)
    def y(r): return r * m.sin(p.mesh_phi/2.0)

    d.vertices.set(    i, [  x(v2D[0]),  y(v2D[0]),  v2D[1]])
    d.vertices.set(k + i, [  x(v2D[0]), -y(v2D[0]),  v2D[1]])

for i in range(n): vset(n, i, v[i])

# Blocks

def bset(i, n, vi, **kwargs):

    k = n + 1

    d.blocks.set(i, vi + [k + i for i in vi], **kwargs)

bset(  0, n, [  0,  1,  2,  3], zone="solid")
bset(  1, n, [  1,  4,  5,  2], zone="solid")
bset(  2, n, [  3,  2,  5,  6], zone="solid")
bset(  3, n, [  4,  7,  8,  9], zone="fluid")
bset(  4, n, [  5,  4,  9, 10], zone="fluid")
bset(  5, n, [  6,  5, 10, 11], zone="fluid")
bset(  6, n, [ 12, 13,  1,  0], zone="vessel")
bset(  7, n, [ 13, 14,  4,  1], zone="vessel")
bset(  8, n, [ 14, 15,  7,  4], zone="vessel")
bset(  9, n, [  7, 15, 16,  8], zone="vessel")
bset( 10, n, [ 17, 18, 13, 12], zone="heater")
bset( 11, n, [ 18, 19, 14, 13], zone="heater")
bset( 12, n, [ 19, 20, 15, 14], zone="heater")
bset( 13, n, [ 20, 21, 22, 15], zone="heater")
bset( 14, n, [ 15, 22, 23, 16], zone="heater")
bset( 15, n, [ 16, 23, 24, 25], zone="heater") # Heater or Space?
bset( 16, n, [  8, 16, 25, 26], zone="space")
bset( 17, n, [  9,  8, 26, 27], zone="space")
bset( 18, n, [ 10,  9, 27, 28], zone="space")
bset( 19, n, [ 11, 10, 28, 29], zone="space")
bset( 20, n, [ 30, 31, 18, 17], zone="space")
bset( 21, n, [ 31, 32, 19, 18], zone="space")
bset( 22, n, [ 32, 33, 20, 19], zone="space")
bset( 23, n, [ 33, 34, 21, 20], zone="space")
bset( 24, n, [ 34, 35, 36, 21], zone="space")
bset( 25, n, [ 21, 36, 37, 22], zone="space")
bset( 26, n, [ 22, 37, 38, 23], zone="space")
bset( 27, n, [ 23, 38, 39, 24], zone="space")
bset( 28, n, [ 24, 39, 40, 41], zone="space")
bset( 29, n, [ 25, 24, 41, 42], zone="space")
bset( 30, n, [ 26, 25, 42, 43], zone="space")
bset( 31, n, [ 27, 26, 43, 44], zone="space")
bset( 32, n, [ 28, 27, 44, 45], zone="space")
bset( 33, n, [ 29, 28, 45, 46], zone="space")

# --------------------------------------------------------------------------- #
# --- blockMeshDict --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

d.header(p.geo_scale)

# --------------------------------------------------------------------------- #

if d.subDict("vertices"):

    d.vertices.write()

if d.subDict("blocks"):

    d.blocks.write()

if d.subDict("edges"):

    pass

if d.subDict("boundary"):

    pass

if d.subDict("mergePatchPairs"):

    pass

# --------------------------------------------------------------------------- #

d.footer()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
