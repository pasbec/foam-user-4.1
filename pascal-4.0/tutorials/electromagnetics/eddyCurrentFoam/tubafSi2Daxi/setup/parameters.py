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

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

geo_scale  = 1e-3

geo_R0 =  53.0
geo_R1 =  62.0
geo_R2 =  70.0

geo_Z0 = -34.0
geo_Z1 =  -8.0
geo_Z2 =   0.0
geo_Z3 =  22.0
geo_Z4 =  27.095
geo_Z5 =  30.0
geo_Z6 =  70.0
geo_Z7 = 142.0

geo_alpha = m.atan(geo_Z5/geo_R0)
geo_Alpha = geo_alpha * 180.0/m.pi

geo_Z     = geo_Z3                      # Solid height
geo_R     = geo_Z3/m.tan(geo_alpha)     # Solid wall contact radius
geo_S     = m.sqrt(geo_Z**2 + geo_R**2) # Solid wall length

# --------------------------------------------------------------------------- #

def csr(): return np.array([1.0, 0.0])
def csz(): return np.array([0.0, 1.0])
def csa(a): return np.array([m.cos(a), m.sin(a)])

mesh_RA  = 0.01

mesh_O0  = mesh_RA * csr()

mesh_Z1  = geo_Z * csz()
mesh_Z0  = 0.5 * mesh_Z1

mesh_S1  = geo_S * csa(geo_alpha)
mesh_S0  = 0.5 * mesh_S1

mesh_A0  = 1.0/3.0*(mesh_O0 + mesh_Z1 + mesh_S1)
mesh_A1  = 1.0/2.0*(mesh_Z1 + mesh_S1)

# --------------------------------------------------------------------------- #
# --- Directories ----------------------------------------------------------- #
# --------------------------------------------------------------------------- #

dir_case = os.path.realpath(csd + '/' + '..')
dir_constant = os.path.realpath(dir_case + '/' + 'constant')
dir_system = os.path.realpath(dir_case+'/' + 'system')

dir_polyMesh = os.path.realpath(dir_constant+'/' + 'polyMesh')
dir_featureEdgeMesh = os.path.realpath(dir_constant+'/' + 'featureEdgeMesh')
dir_triSurface = os.path.realpath(dir_constant+'/' + 'triSurface')

for d in [dir_polyMesh, dir_triSurface, dir_featureEdgeMesh]:
    if not os.path.exists(d): os.makedirs(d)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
