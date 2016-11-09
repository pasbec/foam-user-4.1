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

from foamTools.ioInfo import objectIndent, objectHeader, objectFooter

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

import parameters as par

coil_alpha = (2.0*m.pi)/par.coil_n

coil_arc_alpha = (m.pi/2.0)/par.coil_arc_n

coil_bun_xy = par.coil_scale*par.coil_dxy / (par.coil_bun_n-1)
coil_bun_z  = par.coil_scale*par.coil_dz  / (par.coil_bun_n-1)

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Define short indented line with line break
def i(iL,cS,eS='\n'): return objectIndent(cS + eS,iLevel=iL)




def transform_rmf(x, y, z, n):

    # Switch coordinates to flip coil into upright position
    x1 = z
    y1 = x
    z1 = y + par.coil_scale*par.coil_z

    # Rotate coil into rmf direction
    x2 = m.cos(n*coil_alpha)*x1 - m.sin(n*coil_alpha)*y1
    y2 = m.sin(n*coil_alpha)*x1 + m.cos(n*coil_alpha)*y1
    z2 = z1

    # Shift coil into rmf position
    xr = x2 + par.coil_scale*par.coil_r * m.cos(n*coil_alpha)
    yr = y2 + par.coil_scale*par.coil_r * m.sin(n*coil_alpha)
    zr = z2

    return xr, yr, zr



def a(n): return n*(par.coil_bun_n-1)



def b(n): return a(n)+1



def transform_bundle(n):

    if (n < b(1)):

        bun_x = par.coil_scale*par.coil_xi + n*coil_bun_xy
        bun_y = par.coil_scale*par.coil_yi + n*coil_bun_xy
        bun_r = par.coil_scale*par.coil_ri + n*coil_bun_xy
        bun_z = 0.0

    elif (n > a(1)) and (n < b(2)):

        bun_x = par.coil_scale*par.coil_xi + a(1)*coil_bun_xy
        bun_y = par.coil_scale*par.coil_yi + a(1)*coil_bun_xy
        bun_r = par.coil_scale*par.coil_ri + a(1)*coil_bun_xy
        bun_z = n*coil_bun_z - a(1)*coil_bun_z

    elif (n > a(2)) and (n < b(3)):

        bun_x = par.coil_scale*par.coil_xi + a(3)*coil_bun_xy - n*coil_bun_xy
        bun_y = par.coil_scale*par.coil_yi + a(3)*coil_bun_xy - n*coil_bun_xy
        bun_r = par.coil_scale*par.coil_ri + a(3)*coil_bun_xy - n*coil_bun_xy
        bun_z = a(1)*coil_bun_z

    elif (n > a(3)) and (n < a(4)):

        bun_x = par.coil_scale*par.coil_xi
        bun_y = par.coil_scale*par.coil_yi
        bun_r = par.coil_scale*par.coil_ri
        bun_z = n*coil_bun_z - a(3)*coil_bun_z

    return bun_x, bun_y, bun_r, bun_z

# --------------------------------------------------------------------------- #
# --- Main program sequence ------------------------------------------------- #
# --------------------------------------------------------------------------- #

for geo in range(par.coil_n):

    # Init points
    points = []

    # Init Points/Edges strings
    p=''
    e=''

    for bundle in range(a(4)):

        bun_x, bun_y, bun_r, bun_z = transform_bundle(bundle)

        # Corner ++
        coil_arc_alpha0 = 0.0
        for arc in range(par.coil_arc_n+1):

            x =  bun_x - bun_r + bun_r * m.cos(coil_arc_alpha0 + arc*coil_arc_alpha)
            y =  bun_y - bun_r + bun_r * m.sin(coil_arc_alpha0 + arc*coil_arc_alpha)
            z =  bun_z

            x, y, z = transform_rmf(x, y, z, geo)

            points.append([x, y, z])

        # Corner -+
        coil_arc_alpha0 = m.pi/2.0
        for arc in range(par.coil_arc_n+1):

            x = -bun_x + bun_r + bun_r * m.cos(coil_arc_alpha0 + arc*coil_arc_alpha)
            y =  bun_y - bun_r + bun_r * m.sin(coil_arc_alpha0 + arc*coil_arc_alpha)
            z =  bun_z

            x, y, z = transform_rmf(x, y, z, geo)

            points.append([x, y, z])

        # Corner --
        coil_arc_alpha0 = m.pi
        for arc in range(par.coil_arc_n+1):

            x = -bun_x + bun_r + bun_r * m.cos(coil_arc_alpha0 + arc*coil_arc_alpha)
            y = -bun_y + bun_r + bun_r * m.sin(coil_arc_alpha0 + arc*coil_arc_alpha)
            z =  bun_z

            x, y, z = transform_rmf(x, y, z, geo)

            points.append([x, y, z])

        # Corner +-
        coil_arc_alpha0 = m.pi/2.0*3.0
        for arc in range(par.coil_arc_n+1):

            x =  bun_x - bun_r + bun_r * m.cos(coil_arc_alpha0 + arc*coil_arc_alpha)
            y = -bun_y + bun_r + bun_r * m.sin(coil_arc_alpha0 + arc*coil_arc_alpha)
            z =  bun_z

            x, y, z = transform_rmf(x, y, z, geo)

            points.append([x, y, z])

        shift = bundle*4*(par.coil_arc_n+1)

        # Points string
        for pointi in range(len(points)-shift):

            pointI = pointi + shift

            point = points[pointI]

            p += i(1, '('+str(point[0])+' '+str(point[1])+' '+str(point[2])+')')

        # Edges string
        for pointi in range(len(points)-1-shift):

            pointI = pointi + shift

            e += i(1, '('+str(pointI)+' '+str(pointI+1)+')')

        e += i(1, '('+str(len(points)-1)+' '+str(shift)+')')

    # Write
    name = csn + str(geo) + '.eMesh'
    path = par.dir_featureEdgeMesh + '/' + name
    with open(path, 'w') as f:

        f.write(objectHeader(name, 'featureEdgeMesh'))

        f.write(i(0, '// points:\n'))
        f.write(i(0, str(len(points))))
        f.write(i(0, '('))
        f.write(p)
        f.write(i(0, ')\n'))

        f.write(i(0, '// edges:\n'))
        f.write(i(0, str(len(points))))
        f.write(i(0, '('))
        f.write(e)
        f.write(i(0, ')\n'))

        f.write(objectFooter())

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
