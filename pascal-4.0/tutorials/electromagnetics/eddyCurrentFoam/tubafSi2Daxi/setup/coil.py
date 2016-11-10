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

from foamTools.coil import path
from foamTools.coil import bundleN, bundle, edgeLoopFromPoints
from foamTools.coil import writeCoilFeatureEdgeMeshes
from foamTools.coil import writeEdgeBiotSavartProperties, writeFrequency

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

import parameters as par

# --------------------------------------------------------------------------- #
# --- Coil creation --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

#print(path['loop'](par.coil_path, scale=par.coil_scale))
#exit()

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

    for bundle_i in range(bundleN[b_shape](b)):

        rb, zb = bundle[b_shape](b, bundle_i, r0, z0, par.coil_scale)

        for path_i in range(par.coil_arc_n):

            x =  rb * m.cos(path_i*1.0/par.coil_arc_n * 2.0*m.pi)
            y =  rb * m.sin(path_i*1.0/par.coil_arc_n * 2.0*m.pi)
            z =  zb

            p.append([x, y, z])

        edges[n] += edgeLoopFromPoints(p, len(points[n]))
        points[n] += p

writeCoilFeatureEdgeMeshes(par.dir_case, names, points, edges)

writeEdgeBiotSavartProperties(par.dir_case, names, par.coil_bundle, par.coil_biotSavart)

writeFrequency(par.dir_case, par.coil_frequency)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
