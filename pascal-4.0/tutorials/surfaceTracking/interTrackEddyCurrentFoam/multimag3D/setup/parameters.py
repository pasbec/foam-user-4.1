#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
# Pascal Beckstein (p.beckstein@hzdr.de)

from __future__ import nested_scopes, generators, division, absolute_import
from __future__ import with_statement, print_function, unicode_literals

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

__name__
__path__ = os.path.realpath(__file__)
__base__ = os.path.basename(__path__)
__dir__ = os.path.dirname(__path__)
__head__ = os.path.splitext(__base__)[0]

sys.path.append(os.environ["FOAM_USER_TOOLS"] + "/" + "python")

import math as m
import numpy as np

# --------------------------------------------------------------------------- #
# --- Geometry -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

geo_scale  = 1e-3
geo_space  = 1.0

geo_r0     =   15.0
geo_r1     =   30.0
geo_r2     =   30.0 * 4.0 * geo_space

geo_z0     =  -30.0 * 4.0 * geo_space + 30.0
geo_z1     =    0.0
geo_z2     =   30.0
geo_z3     =   60.0
geo_z4     =   30.0 * 4.0 * geo_space + 30.0

# --------------------------------------------------------------------------- #
# --- Mesh ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

mesh_normal    = -1

mesh_f     =    1.2
mesh_scale =    1.0 * 2e-1

mesh           = {"normal": mesh_normal}

# --------------------------------------------------------------------------- #
# --- Coils ----------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

coil_scale      = 1e-3

coil_bundle     = {"shape": "rectangle",
                   "n":     10,
                   "r":     45.0,
                   "z":     60.0}

coil_path       = {"shape": "racetrack",
                   "n":     9,
                   "r":     10.0 + coil_bundle["r"]/2.0,
                   "x":     (100.0 + coil_bundle["r"])/2.0,
                   "y":     (350.0 + coil_bundle["r"])/2.0}

coils_n         = 6
coils_step      = 285.0
coils_origin    = [0.0, 0.0, geo_z2]

coils_current   = 450.0
coils_nNonOrto  = 10
coils_frequency = 50.0

# We need a current amplitude of
#
# I = 450 A
#
# for achieving
#
# B = sqrt(Bx^2 + By^2 + Bz^2) = 0.4216 mT
#
# at centroid of cylinder for
#
# Ta = 10^5
#
# with Ta = sigma * omega * B^2 * (H/2)^4 / (2 * rho * nu^2)
#

# --------------------------------------------------------------------------- #
# --- Directories ----------------------------------------------------------- #
# --------------------------------------------------------------------------- #

dir_case = os.path.realpath(__dir__ + "/" + "..")

dir_setup = os.path.realpath(dir_case + "/" + "setup")
dir_results = os.path.realpath(dir_case + "/" + "results")

dir_0 = os.path.realpath(dir_case + "/" + "0")
dir_constant = os.path.realpath(dir_case + "/" + "constant")
dir_system = os.path.realpath(dir_case + "/" + "system")

dir_polyMesh = os.path.realpath(dir_constant + "/" + "polyMesh")
dir_featureEdgeMesh = os.path.realpath(dir_constant + "/" + "featureEdgeMesh")
dir_triSurface = os.path.realpath(dir_constant + "/" + "triSurface")

for d in [dir_polyMesh, dir_triSurface, dir_featureEdgeMesh]:
    if not os.path.exists(d): os.makedirs(d)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
