#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

__path__ = os.path.realpath(__file__)
__dir__ = os.path.dirname(__path__)

sys.path.append(os.environ["FOAM_USER_TOOLS"] + "/" + "python")

import math as m
import numpy as np

from foamTools.coil import (inductorCoils, writeCoilFeatureEdgeMeshes,
                            writeEdgeBiotSavartProperties, writeFrequency)

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

import parameters as par

# --------------------------------------------------------------------------- #
# --- Coil creation --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

coils = inductorCoils("RMF", __name__, par.coil_bundle, par.coil_path,
                      par.coils_current, par.coils_n, par.coils_step,
                      origin=par.coils_origin, axis=2, scale=par.coil_scale)

writeCoilFeatureEdgeMeshes(par.dir_case, coils)
writeEdgeBiotSavartProperties(par.dir_case, coils, par.coils_nNonOrto)
writeFrequency(par.dir_case, par.coils_frequency)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
