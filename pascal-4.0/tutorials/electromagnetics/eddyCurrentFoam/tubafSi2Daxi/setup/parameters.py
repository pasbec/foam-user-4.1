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

geo_alpha = m.atan(geo_Z5/geo_R0)      # Cone angle

geo_Zs = geo_Z3                        # Solid height
geo_Rs = geo_Z3/m.tan(geo_alpha)       # Solid wall contact radius
geo_Ss = m.sqrt(geo_Zs**2 + geo_Rs**2) # Solid wall length

# --------------------------------------------------------------------------- #

mesh_scale = 1.0
mesh_space = 1.0

mesh_Ra    = 0.01 # Axis patch radius

mesh_phi   = 5.0  # Wedge angle

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
