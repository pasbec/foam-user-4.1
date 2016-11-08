#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
# Pascal Beckstein (par.beckstein@hzdr.de)

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

import parameters as par

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

d = blockMeshDict(par.dir_polyMesh + '/' + 'blockMeshDict')

# Vertices

def vmg(v): return np.linalg.norm(v)
def vr0(r): return np.array([r, 0.0])
def v0z(z): return np.array([par.mesh_R['axis'], z])
def vrz(r, z): return np.array([r, z])
def vsa(s, a): return s * np.array([m.cos(a), m.sin(a)])
def vwg(al, wl=None):
    if wl is None: wl = np.ones(len(al))
    a = np.zeros(2)
    for i, ai in enumerate(al): a = np.add(a, wl[i]*ai)
    return a/np.sum(wl)

v      = dict()

v[  0] = v0z(0.0)
v[  6] = v0z(par.geo_Z['solid'])
v[  3] = vwg([v[0], v[6]])                    # Centroid of solid axis height
v[  4] = vrz(par.geo_R['solid'], v[6][1])
v[  1] = vwg([v[0], v[4]])                    # Centroid of solid wall contact
v[  2] = vwg([v[0], v[4], v[6]], [1, 1.3, 1]) # Centroid of solid (weighted)
v[  5] = vwg([v[4], v[6]])                    # Centroid of solid radius width

v[  7] = vrz(par.geo_R[0], par.geo_Z[5])
v[  8] = vrz(v[7][0], par.geo_Z[6])
v[  9] = vrz(v[4][0], v[8][1])
v[ 10] = vrz(v[5][0], v[8][1])
v[ 11] = v0z(v[8][1])

v[ 12] = v0z(par.geo_Z[1])
v[ 13] = v[1] + par.geo_G
v[ 14] = v[4] + par.geo_G
v[ 15] = vrz(par.geo_R[1], par.geo_Z[4])
v[ 16] = vrz(v[15][0], v[8][1])

v[ 17] = v0z(par.geo_Z[0])
v[ 18] = vrz(v[13][0], v[17][1])
v[ 19] = vrz(v[14][0], v[17][1])
v[ 20] = vrz(v[15][0], v[17][1])
v[ 21] = vrz(par.geo_R[2], v[20][1])
v[ 22] = vrz(v[21][0], v[15][1])
v[ 23] = vrz(v[21][0], v[8][1])
v[ 24] = vrz(v[21][0], par.geo_Z[7])
v[ 25] = vrz(v[15][0], v[24][1])
v[ 26] = vrz(v[7][0], v[24][1])
v[ 27] = vrz(v[9][0], v[24][1])
v[ 28] = vrz(v[5][0], v[24][1])
v[ 29] = v0z(v[28][1])

v[ 30] = v0z(par.mesh_Z['C'] - par.mesh_Z['inf'])
v[ 31] = vrz(v[18][0], v[30][1])
v[ 32] = vrz(v[19][0], v[30][1])
v[ 33] = vrz(v[20][0], v[30][1])
v[ 34] = vrz(v[21][0], v[30][1])
v[ 35] = vrz(par.mesh_R['inf'], v[30][1])
v[ 36] = vrz(v[35][0], v[21][1])
v[ 37] = vrz(v[36][0], v[22][1])
v[ 38] = vrz(v[36][0], v[23][1])
v[ 39] = vrz(v[36][0], v[24][1])
v[ 40] = vrz(v[36][0], par.mesh_Z['C'] + par.mesh_Z['inf'])
v[ 41] = vrz(v[24][0], v[40][1])
v[ 42] = vrz(v[25][0], v[40][1])
v[ 43] = vrz(v[26][0], v[40][1])
v[ 44] = vrz(v[27][0], v[40][1])
v[ 45] = vrz(v[28][0], v[40][1])
v[ 46] = vrz(v[29][0], v[40][1])

nv = len(v)

def vset(nv, i, v2D):

    kv = nv + 1

    def x(r): return r * m.cos(par.mesh_phi/2.0)
    def y(r): return r * m.sin(par.mesh_phi/2.0)

    d.vertices.set(     i, [  x(v2D[0]),  y(v2D[0]),  v2D[1]])
    d.vertices.set(kv + i, [  x(v2D[0]), -y(v2D[0]),  v2D[1]])

for i in range(nv): vset(nv, i, v[i])

# Blocks

def bset(i, nv, vi, **kwargs):

    kv = nv + 1

    d.blocks.set(i, vi + [kv + i for i in vi], **kwargs)

bset(  0, nv, [  0,  1,  2,  3], zone="solid")
bset(  1, nv, [  1,  4,  5,  2], zone="solid")
bset(  2, nv, [  3,  2,  5,  6], zone="solid")
bset(  3, nv, [  4,  7,  8,  9], zone="fluid")
bset(  4, nv, [  5,  4,  9, 10], zone="fluid")
bset(  5, nv, [  6,  5, 10, 11], zone="fluid")
bset(  6, nv, [ 12, 13,  1,  0], zone="vessel")
bset(  7, nv, [ 13, 14,  4,  1], zone="vessel")
bset(  8, nv, [ 14, 15,  7,  4], zone="vessel")
bset(  9, nv, [  7, 15, 16,  8], zone="vessel")
bset( 10, nv, [ 17, 18, 13, 12], zone="heater")
bset( 11, nv, [ 18, 19, 14, 13], zone="heater")
bset( 12, nv, [ 19, 20, 15, 14], zone="heater")
bset( 13, nv, [ 20, 21, 22, 15], zone="heater")
bset( 14, nv, [ 15, 22, 23, 16], zone="heater")
bset( 15, nv, [ 16, 23, 24, 25], zone="heater") # Heater or Space?
bset( 16, nv, [  8, 16, 25, 26], zone="space")
bset( 17, nv, [  9,  8, 26, 27], zone="space")
bset( 18, nv, [ 10,  9, 27, 28], zone="space")
bset( 19, nv, [ 11, 10, 28, 29], zone="space")
bset( 20, nv, [ 30, 31, 18, 17], zone="space")
bset( 21, nv, [ 31, 32, 19, 18], zone="space")
bset( 22, nv, [ 32, 33, 20, 19], zone="space")
bset( 23, nv, [ 33, 34, 21, 20], zone="space")
bset( 24, nv, [ 34, 35, 36, 21], zone="space")
bset( 25, nv, [ 21, 36, 37, 22], zone="space")
bset( 26, nv, [ 22, 37, 38, 23], zone="space")
bset( 27, nv, [ 23, 38, 39, 24], zone="space")
bset( 28, nv, [ 24, 39, 40, 41], zone="space")
bset( 29, nv, [ 25, 24, 41, 42], zone="space")
bset( 30, nv, [ 26, 25, 42, 43], zone="space")
bset( 31, nv, [ 27, 26, 43, 44], zone="space")
bset( 32, nv, [ 28, 27, 44, 45], zone="space")
bset( 33, nv, [ 29, 28, 45, 46], zone="space")

# Distributions

def n(l): s = 0.1; return int(m.ceil(abs(s*par.mesh_scale*l)))

nr                 = dict()
nr['space_side']   = n(v[35][0] - v[34][0])
nr['gap']          = n(v[42][0] - v[43][0])
nr['heater_side']  = n(v[41][0] - v[42][0])
nr['heater_width'] = n(v[19][0] - v[18][0])
nr['fluid_width0'] = n(v[10][0] - v[11][0])
nr['fluid_width1'] = n(v[9][0] - v[10][0])
nr['fluid_width2'] = n(v[8][0] - v[9][0])

print(nr)

nz                   = dict()
nz['space_top']      = n(v[40][1] - v[39][1])
nz['space_bottom']   = n(v[36][1] - v[35][1])
nz['space_internal'] = n(v[29][1] - v[11][1])
nz['heater_bottom']  = n(v[37][1] - v[36][1])
nz['fluid_height']   = n(v[11][1] - v[6][1])

print(nz)

d.blocks.distribution.set( 24, "x", nr['space_side'])
d.blocks.distribution.set( 30, "x", nr['gap'])
d.blocks.distribution.set( 29, "x", nr['heater_side'])
d.blocks.distribution.set( 11, "x", nr['heater_width'])
d.blocks.distribution.set(  5, "x", nr['fluid_width0'])
d.blocks.distribution.set(  4, "x", nr['fluid_width1'])
d.blocks.distribution.set(  3, "x", nr['fluid_width2'])

d.blocks.distribution.set( 33, "y", nz['space_top'])
d.blocks.distribution.set( 20, "y", nz['space_bottom'])
d.blocks.distribution.set( 19, "y", nz['space_internal'])
d.blocks.distribution.set( 10, "y", nz['heater_bottom'])
d.blocks.distribution.set(  5, "y", nz['fluid_height'])

# Boundary faces

d.boundaryFaces.set(  0, "axis", [  0,  2,  5,  6, 10, 19, 20, 33], "x-")
d.boundaryFaces.set(  1, "front", [ i for i in range(len(d.blocks.labels)) ], "z-")
d.boundaryFaces.set(  2, "back", [ i for i in range(len(d.blocks.labels)) ], "z+")
d.boundaryFaces.set(  3, "infinity", [ 24, 25, 26, 27, 28], "x+")
d.boundaryFaces.set(  4, "infinity", [ 28, 29, 30, 31, 32, 33], "y+")
d.boundaryFaces.set(  5, "infinity", [ 20, 21, 22, 23, 24], "y-")

# --------------------------------------------------------------------------- #
# --- blockMeshDict --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

d.header(par.geo_scale)

# --------------------------------------------------------------------------- #

if d.subDict("vertices"):

    d.vertices.write()

if d.subDict("blocks"):

    d.blocks.write()

if d.subDict("edges"):

    pass

if d.subDict("boundary"):

    if d.boundarySubDict("axis", "patch"):

        d.boundaryFaces.write()

    if d.boundarySubDict("front", "wedge"):

        d.boundaryFaces.write()

    if d.boundarySubDict("back", "wedge"):

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
