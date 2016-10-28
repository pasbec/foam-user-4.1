#!/usr/bin/python
# July 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

from foamTools.ioInfo import fileGetPath, objectIndent, objectHeader, objectFooter

import math as m

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

geo_x = 0.065
geo_y = 0.120
geo_r = 0.008
geo_z = 0.009

off_n = 3
off_z = 0.013

arc_n = 12
arc_alpha = (m.pi/2.0)/arc_n

fileNameBase = 'coil'

# --------------------------------------------------------------------------- #

# Define short indented line with line break
def i(iL,cS,eS='\n'): return objectIndent(cS + eS,iLevel=iL)

# --------------------------------------------------------------------------- #

for o in range(off_n):

    # Init points
    points = []

    # Init Points/Edges strings
    p=''
    e=''

    # Corner ++
    arc_alpha0 = 0.0
    for n in range(arc_n+1):

        x =  geo_x - geo_r + geo_r * m.cos(arc_alpha0 + n*arc_alpha)
        y =  geo_y - geo_r + geo_r * m.sin(arc_alpha0 + n*arc_alpha)
        z =  geo_z + o*off_z

        points.append([x, y, z])

    # Corner -+
    arc_alpha0 = m.pi/2.0
    for n in range(arc_n+1):

        x = -geo_x + geo_r + geo_r * m.cos(arc_alpha0 + n*arc_alpha)
        y =  geo_y - geo_r + geo_r * m.sin(arc_alpha0 + n*arc_alpha)
        z =  geo_z + o*off_z

        points.append([x, y, z])

    # Corner --
    arc_alpha0 = m.pi
    for n in range(arc_n+1):

        x = -geo_x + geo_r + geo_r * m.cos(arc_alpha0 + n*arc_alpha)
        y = -geo_y + geo_r + geo_r * m.sin(arc_alpha0 + n*arc_alpha)
        z =  geo_z + o*off_z

        points.append([x, y, z])

    # Corner +-
    arc_alpha0 = m.pi/2.0*3.0
    for n in range(arc_n+1):

        x =  geo_x - geo_r + geo_r * m.cos(arc_alpha0 + n*arc_alpha)
        y = -geo_y + geo_r + geo_r * m.sin(arc_alpha0 + n*arc_alpha)
        z =  geo_z + o*off_z

        points.append([x, y, z])

    # Points string
    for pointI in range(len(points)):

        point = points[pointI]

        p += i(1, '('+str(point[0])+' '+str(point[1])+' '+str(point[2])+')')

    # Edges string
    for pointI in range(len(points)-1):

        e += i(1, '('+str(pointI)+' '+str(pointI+1)+')')

    e += i(1, '('+str(len(points)-1)+' '+str(0)+')')

    # Write
    fileName = fileNameBase + str(o) + '.eMesh'
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
