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

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

coil_scale = 1

coil_xi    = 0.050 # half inner width:  100/2.0
coil_yi    = 0.175 # half inner height: 350/2.0
#coil_y     = 500.0 # half inner height: 1000000/2.0
coil_ri    = 0.010 # inner radius
coil_dxy   = 0.045 # thickness
coil_dz    = 0.060 # height

coil_r     = 0.285 - coil_dz/2.0
coil_z     = 0.030

coil_n     = 6
coil_arc_n = 9
coil_bun_n = 10 # 4*10 - 2 = 38

# --------------------------------------------------------------------------- #

geo_scale  = 1e-3

geo_r0     =   15.0
geo_r1     =   30.0
#geo_r2     =  120.0
geo_r2     =  500.0

#geo_z0     =  -60.0
geo_z0     = -440.0
geo_z1     =    0.0
geo_z2     =   30.0
geo_z3     =   60.0
#geo_z4     =  120.0
geo_z4     =  500.0

# --------------------------------------------------------------------------- #

mesh_f     = 1.2
mesh_scale = 0.5

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
