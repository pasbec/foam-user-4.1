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
v[  4] = vrz(p.geo_R['solid'], p.geo_Z['solid'])
v[  1] = vct(v[0], v[4])
v[  2] = vct(v[0], v[4], v[6])
v[  5] = vct(v[4], v[6])

v[  7] = vrz(p.geo_R[0], p.geo_Z[5])
v[  8] = vrz(p.geo_R[0], p.geo_Z[6])
v[  9] = vrz(v[4][0], p.geo_Z[6])
v[ 10] = vrz(v[5][0], p.geo_Z[6])
v[ 11] = v0z(p.geo_Z[6])

v[ 12] = v0z(p.geo_Z[1])
v[ 13] = v[1] + p.geo_G
v[ 14] = v[4] + p.geo_G
v[ 15] = vrz(p.geo_R[1], p.geo_Z[4])
v[ 16] = vrz(p.geo_R[1], p.geo_Z[6])

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
bset(  6, n, [ 12, 13,  1,  0], zone="container")
bset(  7, n, [ 13, 14,  4,  1], zone="container")
bset(  8, n, [ 14, 15,  7,  4], zone="container")

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
