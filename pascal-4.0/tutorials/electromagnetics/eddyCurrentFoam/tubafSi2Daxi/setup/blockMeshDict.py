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

from foamTools.expansion import expansion_de_e, expansion_n_ds
from foamTools.blockMeshDict import blockMeshDict

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

import parameters as p

axi_Phi        = 5.0
axi_phi        = m.pi * axi_Phi/180.0
axi_phiHalf    = axi_phi/2.0
axi_phiHalfCos = m.cos(axi_phiHalf)
axi_phiHalfSin = m.sin(axi_phiHalf)

def x(r): return r * axi_phiHalfCos
def y(r): return r * axi_phiHalfSin
def ap(v2D): return [  x(v2D[0]),  y(v2D[0]),  v2D[1]]
def an(v2D): return [  x(v2D[0]), -y(v2D[0]),  v2D[1]]

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

d = blockMeshDict(p.dir_polyMesh + '/' + 'blockMeshDict')

# Vertices

d.vertices.set(  0, ap(p.mesh_O0))
d.vertices.set(  1, ap(p.mesh_S0))
d.vertices.set(  2, ap(p.mesh_A0))
d.vertices.set(  3, ap(p.mesh_Z0))
d.vertices.set(  4, ap(p.mesh_S1))
d.vertices.set(  5, ap(p.mesh_A1))
d.vertices.set(  6, ap(p.mesh_Z1))

d.vertices.set(100, an(p.mesh_O0))
d.vertices.set(101, an(p.mesh_S0))
d.vertices.set(102, an(p.mesh_A0))
d.vertices.set(103, an(p.mesh_Z0))
d.vertices.set(104, an(p.mesh_S1))
d.vertices.set(105, an(p.mesh_A1))
d.vertices.set(106, an(p.mesh_Z1))

# Blocks

d.blocks.set(  0, [   0,   1,   2,   3, 100, 101, 102, 103], zone="solid")
d.blocks.set(  1, [   1,   4,   5,   2, 101, 104, 105, 102], zone="solid")
d.blocks.set(  2, [   3,   2,   5,   6, 103, 102, 105, 106], zone="solid")

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
