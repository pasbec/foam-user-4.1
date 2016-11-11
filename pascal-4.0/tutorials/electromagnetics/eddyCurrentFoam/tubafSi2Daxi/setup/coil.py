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

from foamTools.coil import coil
from foamTools.coil import writeCoilFeatureEdgeMeshes
from foamTools.coil import writeEdgeBiotSavartProperties, writeFrequency

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

import parameters as par

# --------------------------------------------------------------------------- #
# --- Coil creation --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

coils = dict()

for n in range(par.coil_n):

    coils[n] = coil(csn + str(n),
                    par.coil_reverse, par.coil_current, par.coil_phase,
                    par.coil_bundle, par.coil_path,
                    translate=[0.0, 0.0, par.coil_z + n*par.coil_dz],
                    scale=par.coil_scale)

writeCoilFeatureEdgeMeshes(par.dir_case, coils)
writeEdgeBiotSavartProperties(par.dir_case, coils, par.coil_nNonOrto)
writeFrequency(par.dir_case, par.coil_frequency)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
