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

    r0 = par.coil_r
    z0 = par.coil_z + n*par.coil_dn

    b       = par.coil_bundle
    b_shape = b['shape']

    for bun_i in range(bundleN[b_shape](b)):

        p = list()

        rb, zb = bundle[b_shape](b, bun_i, r0, z0, par.coil_scale)

        for arc_i in range(par.coil_arc_n):

            x =  rb * m.cos(arc_i*2*m.pi/par.coil_arc_n)
            y =  rb * m.sin(arc_i*2*m.pi/par.coil_arc_n)
            z =  zb

            p.append([x, y, z])

        edges[n] += edgeLoopFromPoints(p, len(points[n]))
        points[n] += p

    writeCoilFeatureEdgeMesh(par.dir_featureEdgeMesh,
                             names[n], points[n], edges[n])

writeEdgeBiotSavartProperties(par.dir_constant, names,
                              par.coil_bundle, par.coil_biotSavart)

writeFrequency(par.dir_constant, par.coil_frequency)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
