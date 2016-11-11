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

geo_scale      = 1e-3

geo_R          = dict()
geo_R[0]       =  53.0
geo_R[1]       =  62.0
geo_R[2]       =  70.0

geo_Z          = dict()
geo_Z[0]       = -34.0
geo_Z[1]       =  -8.0
geo_Z[2]       =   0.0
geo_Z[3]       =  22.0
geo_Z[4]       =  27.095
geo_Z[5]       =  30.0
geo_Z[6]       =  70.0
geo_Z[7]       = 142.0

geo_alpha      = m.atan(geo_Z[5]/geo_R[0])                # Cone angle (rad)

geo_magG       = (geo_Z[2] - geo_Z[1]) * m.cos(geo_alpha) # Non-conducting gap size
geo_G          = geo_magG * np.array([m.sin(geo_alpha), -m.cos(geo_alpha)])

geo_Z['solid'] =  geo_Z[3]                                # Solid height
geo_R['solid'] =  geo_Z[3] / m.tan(geo_alpha)             # Solid wall contact radius

# --------------------------------------------------------------------------- #

mesh_scale     = 1.0
mesh_space     = 1.0

mesh_normal    = 1
mesh_thickness = 10.0
mesh_angle     =  5.0

mesh_R         = dict()
mesh_R['axis'] = 1.0                                      # Axis patch radius
mesh_R['inf']  = 4.0 * mesh_space * geo_R[2]              # Infinity patch radius

mesh_Z         = dict()
mesh_Z['C']    = geo_Z[0] + 0.5*(geo_Z[7] - geo_Z[0])
mesh_Z['inf']  = 4.0 * mesh_space * (geo_Z[6] - geo_Z[3]) # Infinity patch radius

mesh_R['inf']  = max(mesh_R['inf'], mesh_Z['inf'])
mesh_Z['inf']  = max(mesh_Z['inf'], mesh_R['inf'])

# --------------------------------------------------------------------------- #

coil_scale      = 1e-3

coil_z          =  7.0
coil_dz         = 14.9

coil_n          = 10

coil_path       = {'shape': 'loop',
                   'n': 36,
                   'r': 105.0}

#coil_bundle     = {'shape': 'point'}
#coil_bundle     = {'shape': 'circle',
                   #'r':     5.0,
                   #'n':     10}
coil_bundle     = {'shape': 'rectangle',
                   'n':     10,
                   'x':     10.0,
                   'y':     8.0}

coil_reverse    = False
coil_current    = m.sqrt(2.0) * 260.0
coil_phase      = 0.0
coil_nNonOrto   = 10
coil_frequency  = 6300.0

# --------------------------------------------------------------------------- #
# --- Directories ----------------------------------------------------------- #
# --------------------------------------------------------------------------- #

dir_case = os.path.realpath(csd + '/' + '..')
dir_constant = os.path.realpath(dir_case + '/' + 'constant')
dir_system = os.path.realpath(dir_case + '/' + 'system')

dir_polyMesh = os.path.realpath(dir_constant + '/' + 'polyMesh')
dir_featureEdgeMesh = os.path.realpath(dir_constant + '/' + 'featureEdgeMesh')
dir_triSurface = os.path.realpath(dir_constant + '/' + 'triSurface')

for d in [dir_polyMesh, dir_triSurface, dir_featureEdgeMesh]:
    if not os.path.exists(d): os.makedirs(d)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
