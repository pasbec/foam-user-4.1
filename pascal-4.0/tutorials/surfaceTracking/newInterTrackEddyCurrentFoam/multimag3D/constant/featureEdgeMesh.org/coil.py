#!/usr/bin/python
# July 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

from foamTools.ioInfo import fileGetPath, objectIndent, objectHeader, objectFooter

import math as m

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

geo_x = 0.080
geo_y = 0.100
geo_r = 0.020
geo_d = 0.020
geo_z = 0.000

rmf_r = 0.300
rmf_z = 0.030
rmf_n = 6
rmf_alpha = (2.0*m.pi)/rmf_n

arc_n = 9
arc_alpha = (m.pi/2.0)/arc_n

bundle_n = 10 # 4*10 - 2 = 38
bundle_r = geo_r/(bundle_n-1)
bundle_z = geo_d/(bundle_n-1)

fileNameBase = 'coil'

# --------------------------------------------------------------------------- #

# Define short indented line with line break
def i(iL,cS,eS='\n'): return objectIndent(cS + eS,iLevel=iL)

# --------------------------------------------------------------------------- #

def transform_rmf(x, y, z, n):

    # Switch coordinates to flip coil into upright position
    x1 = z
    y1 = x
    z1 = y + rmf_z

    # Rotate coil into rmf direction
    x2 = m.cos(n*rmf_alpha)*x1 - m.sin(n*rmf_alpha)*y1
    y2 = m.sin(n*rmf_alpha)*x1 + m.cos(n*rmf_alpha)*y1
    z2 = z1

    # Shift coil into rmf position
    xr = x2 + rmf_r * m.cos(n*rmf_alpha)
    yr = y2 + rmf_r * m.sin(n*rmf_alpha)
    zr = z2

    return xr, yr, zr

for rmf in range(rmf_n):

    # Init points
    points = []

    # Init Points/Edges strings
    p=''
    e=''

    def a(n): return n*(bundle_n-1)
    def b(n): return a(n)+1

    for bundle in range(a(4)):

        if (bundle < b(1)):

            bun_x = geo_x + bundle*bundle_r
            bun_y = geo_y + bundle*bundle_r
            bun_r = geo_r + bundle*bundle_r
            bun_z = geo_z

        elif (bundle > a(1)) and (bundle < b(2)):

            bun_x = geo_x + a(1)*bundle_r
            bun_y = geo_y + a(1)*bundle_r
            bun_r = geo_r + a(1)*bundle_r
            bun_z = geo_z + bundle*bundle_z - a(1)*bundle_z

        elif (bundle > a(2)) and (bundle < b(3)):

            bun_x = geo_x + a(3)*bundle_r - bundle*bundle_r
            bun_y = geo_y + a(3)*bundle_r - bundle*bundle_r
            bun_r = geo_r + a(3)*bundle_r - bundle*bundle_r
            bun_z = geo_z + a(1)*bundle_z

        elif (bundle > a(3)) and (bundle < a(4)):

            bun_x = geo_x
            bun_y = geo_y
            bun_r = geo_r
            bun_z = geo_z + bundle*bundle_z - a(3)*bundle_z

        # Corner ++
        arc_alpha0 = 0.0
        for arc in range(arc_n+1):

            x =  bun_x - bun_r + bun_r * m.cos(arc_alpha0 + arc*arc_alpha)
            y =  bun_y - bun_r + bun_r * m.sin(arc_alpha0 + arc*arc_alpha)
            z =  bun_z

            x, y, z = transform_rmf(x, y, z, rmf)

            points.append([x, y, z])

        # Corner -+
        arc_alpha0 = m.pi/2.0
        for arc in range(arc_n+1):

            x = -bun_x + bun_r + bun_r * m.cos(arc_alpha0 + arc*arc_alpha)
            y =  bun_y - bun_r + bun_r * m.sin(arc_alpha0 + arc*arc_alpha)
            z =  bun_z

            x, y, z = transform_rmf(x, y, z, rmf)

            points.append([x, y, z])

        # Corner --
        arc_alpha0 = m.pi
        for arc in range(arc_n+1):

            x = -bun_x + bun_r + bun_r * m.cos(arc_alpha0 + arc*arc_alpha)
            y = -bun_y + bun_r + bun_r * m.sin(arc_alpha0 + arc*arc_alpha)
            z =  bun_z

            x, y, z = transform_rmf(x, y, z, rmf)

            points.append([x, y, z])

        # Corner +-
        arc_alpha0 = m.pi/2.0*3.0
        for arc in range(arc_n+1):

            x =  bun_x - bun_r + bun_r * m.cos(arc_alpha0 + arc*arc_alpha)
            y = -bun_y + bun_r + bun_r * m.sin(arc_alpha0 + arc*arc_alpha)
            z =  bun_z

            x, y, z = transform_rmf(x, y, z, rmf)

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
    fileName = fileNameBase + str(rmf) + '.eMesh'
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
