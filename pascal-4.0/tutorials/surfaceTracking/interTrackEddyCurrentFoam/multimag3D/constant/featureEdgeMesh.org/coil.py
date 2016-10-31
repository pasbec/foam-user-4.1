#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import math as m

from foamTools.ioInfo import fileGetPath, objectIndent, objectHeader, objectFooter

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

coil_x   = 0.050 # half inner width:  100/2.0
coil_y   = 0.175 # half inner height: 350/2.0
#coil_y   = 500.0 # half inner height: 1000000/2.0
coil_r   = 0.010 # inner radius
coil_dxy = 0.045 # thickness
coil_dz  = 0.060 # height

geo_r = 0.285 - coil_dz/2.0
geo_z = 0.030
geo_n = 6
geo_alpha = (2.0*m.pi)/geo_n

arc_n = 9
arc_alpha = (m.pi/2.0)/arc_n

bundle_n  = 10 # 4*10 - 2 = 38
bundle_xy = coil_dxy / (bundle_n-1)
bundle_z  = coil_dz  / (bundle_n-1)

fileNameBase = 'coil'

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Define short indented line with line break
def i(iL,cS,eS='\n'): return objectIndent(cS + eS,iLevel=iL)




def transform_rmf(x, y, z, n):

    # Switch coordinates to flip coil into upright position
    x1 = z
    y1 = x
    z1 = y + geo_z

    # Rotate coil into rmf direction
    x2 = m.cos(n*geo_alpha)*x1 - m.sin(n*geo_alpha)*y1
    y2 = m.sin(n*geo_alpha)*x1 + m.cos(n*geo_alpha)*y1
    z2 = z1

    # Shift coil into rmf position
    xr = x2 + geo_r * m.cos(n*geo_alpha)
    yr = y2 + geo_r * m.sin(n*geo_alpha)
    zr = z2

    return xr, yr, zr



def a(n): return n*(bundle_n-1)



def b(n): return a(n)+1



def transform_bundle(n):

    if (n < b(1)):

        bun_x = coil_x + n*bundle_xy
        bun_y = coil_y + n*bundle_xy
        bun_r = coil_r + n*bundle_xy
        bun_z = 0.0

    elif (n > a(1)) and (n < b(2)):

        bun_x = coil_x + a(1)*bundle_xy
        bun_y = coil_y + a(1)*bundle_xy
        bun_r = coil_r + a(1)*bundle_xy
        bun_z = n*bundle_z - a(1)*bundle_z

    elif (n > a(2)) and (n < b(3)):

        bun_x = coil_x + a(3)*bundle_xy - n*bundle_xy
        bun_y = coil_y + a(3)*bundle_xy - n*bundle_xy
        bun_r = coil_r + a(3)*bundle_xy - n*bundle_xy
        bun_z = a(1)*bundle_z

    elif (n > a(3)) and (n < a(4)):

        bun_x = coil_x
        bun_y = coil_y
        bun_r = coil_r
        bun_z = n*bundle_z - a(3)*bundle_z

    return bun_x, bun_y, bun_r, bun_z

# --------------------------------------------------------------------------- #
# --- Main program sequence ------------------------------------------------- #
# --------------------------------------------------------------------------- #

for geo in range(geo_n):

    # Init points
    points = []

    # Init Points/Edges strings
    p=''
    e=''

    for bundle in range(a(4)):

        bun_x, bun_y, bun_r, bun_z = transform_bundle(bundle)

        # Corner ++
        arc_alpha0 = 0.0
        for arc in range(arc_n+1):

            x =  bun_x - bun_r + bun_r * m.cos(arc_alpha0 + arc*arc_alpha)
            y =  bun_y - bun_r + bun_r * m.sin(arc_alpha0 + arc*arc_alpha)
            z =  bun_z

            x, y, z = transform_rmf(x, y, z, geo)

            points.append([x, y, z])

        # Corner -+
        arc_alpha0 = m.pi/2.0
        for arc in range(arc_n+1):

            x = -bun_x + bun_r + bun_r * m.cos(arc_alpha0 + arc*arc_alpha)
            y =  bun_y - bun_r + bun_r * m.sin(arc_alpha0 + arc*arc_alpha)
            z =  bun_z

            x, y, z = transform_rmf(x, y, z, geo)

            points.append([x, y, z])

        # Corner --
        arc_alpha0 = m.pi
        for arc in range(arc_n+1):

            x = -bun_x + bun_r + bun_r * m.cos(arc_alpha0 + arc*arc_alpha)
            y = -bun_y + bun_r + bun_r * m.sin(arc_alpha0 + arc*arc_alpha)
            z =  bun_z

            x, y, z = transform_rmf(x, y, z, geo)

            points.append([x, y, z])

        # Corner +-
        arc_alpha0 = m.pi/2.0*3.0
        for arc in range(arc_n+1):

            x =  bun_x - bun_r + bun_r * m.cos(arc_alpha0 + arc*arc_alpha)
            y = -bun_y + bun_r + bun_r * m.sin(arc_alpha0 + arc*arc_alpha)
            z =  bun_z

            x, y, z = transform_rmf(x, y, z, geo)

            points.append([x, y, z])

        shift = bundle*4*(arc_n+1)

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
    fileName = fileNameBase + str(geo) + '.eMesh'
    with open(fileGetPath(fileName),'w') as f:

        f.write(objectHeader(fileName, 'featureEdgeMesh'))

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
