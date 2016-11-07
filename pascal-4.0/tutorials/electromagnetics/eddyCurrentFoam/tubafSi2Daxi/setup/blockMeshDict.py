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

def cr0(r): return np.array([r, 0.0])
def c0z(z): return np.array([p.mesh_Ra, z])
def crz(r, z): return np.array([r, z])
def cla(l, a): return l * np.array([m.cos(a), m.sin(a)])

v     = dict()
v[  0] = c0z(0.0)
v[  6] = c0z(p.geo_Zs)
v[  3] = 0.5 * v[6]
v[  4] = crz(p.geo_Rs, p.geo_Zs)
v[  1] = 0.5 * v[4]
v[  2] = 1.0/3.0 * (v[0] + v[4] + v[6])
v[  5] = 1.0/2.0 * (v[4] + v[6])
v[  7] = crz(p.geo_R0, p.geo_Z5)
v[  8] = crz(p.geo_R0, p.geo_Z6)
v[  9] = crz(v[4][0], p.geo_Z6)
v[ 10] = crz(v[5][0], p.geo_Z6)
v[ 11] = c0z(p.geo_Z6)

n = len(v)

def vset(n, i, v2D):

    k = n + 1

    def x(r): return r * m.cos(m.pi/180.0 * p.mesh_phi/2.0)
    def y(r): return r * m.sin(m.pi/180.0 * p.mesh_phi/2.0)

    d.vertices.set(    i, [  x(v2D[0]),  y(v2D[0]),  v2D[1]])
    d.vertices.set(k + i, [  x(v2D[0]), -y(v2D[0]),  v2D[1]])

for i in range(n): vset(n, i, v[i])

# Blocks

def bset(i, n, vi, **kwargs):

    k = n + 1

    d.blocks.set(i, vi + [k + i for i in vi], **kwargs)

bset( 0, n, [  0,  1,  2,  3], zone="solid")
bset( 1, n, [  1,  4,  5,  2], zone="solid")
bset( 2, n, [  3,  2,  5,  6], zone="solid")
bset( 3, n, [  4,  7,  8,  9], zone="fluid")
bset( 4, n, [  5,  4,  9, 10], zone="fluid")
bset( 5, n, [  6,  5, 10, 11], zone="fluid")

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
