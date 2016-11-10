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

from foamTools.coil import bundleN, bundle, edgeLoopFromPoints
from foamTools.coil import writeCoilFeatureEdgeMesh
from foamTools.coil import writeEdgeBiotSavartProperties, writeFrequency

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

import parameters as par

coil_r      = par.coil_scale*par.coil_r
coil_z      = par.coil_scale*par.coil_z

coil_dr     = par.coil_scale*par.coil_dr
coil_dz     = par.coil_scale*par.coil_dz

coil_dn     = par.coil_scale*par.coil_dn

# --------------------------------------------------------------------------- #
# --- Coil creation --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

names  = dict()
points = dict()
edges  = dict()

for n in range(par.coil_n):

    names[n]  = csn + str(n)
    points[n] = list()
    edges[n]  = list()

    r0 = coil_r
    z0 = coil_z + n*coil_dn

    for bun_i in range(bundleN['circle'](par.coil_bun_n)):

        p = list()

        rb, zb = bundle['rectangle'](bun_i, par.coil_bun_n, coil_dr, coil_dz, r0, z0)

        for arc_i in range(par.coil_arc_n):

            x =  rb * m.cos(arc_i*2*m.pi/par.coil_arc_n)
            y =  rb * m.sin(arc_i*2*m.pi/par.coil_arc_n)
            z =  zb

            p.append([x, y, z])

        edges[n] += edgeLoopFromPoints(p, len(points[n]))
        points[n] += p

    writeCoilFeatureEdgeMesh(par.dir_featureEdgeMesh,
                             names[n], points[n], edges[n])

writeEdgeBiotSavartProperties(par.dir_constant, par.coil_nNonOrth, names,
                              par.coil_reverse, par.coil_current, par.coil_phase)

writeFrequency(par.dir_constant, par.coil_frequency)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
