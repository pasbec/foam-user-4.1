#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Module template
# March 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

from __future__ import nested_scopes, generators, division, absolute_import
from __future__ import with_statement, print_function, unicode_literals

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

__name__
__path__ = os.path.realpath(__file__)
__base__ = os.path.basename(__path__)
__dir__ = os.path.dirname(__path__)
__head__ = os.path.splitext(__base__)[0]

sys.path.append("/usr/lib/freecad/lib")

import math as m
import numpy as np

import FreeCAD, Sketcher, Draft, Part, PartDesign, Mesh, MeshPart

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# --- Functions ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def addPolyLine(s, v, v2D):

    l = len(v) - 1

    def V(i): return FreeCAD.Vector(v2D[i][0], v2D[i][1], 0.0)

    for i in range(l):

        s.addGeometry(Part.Line(V(v[i]), V(v[i+1])))

    s.addGeometry(Part.Line(V(v[l]), V(v[0])))



def faceShell(shd):

    dList = list([shd]) if isinstance(shd, tuple) else shd

    if not isinstance(dList, list):

        raise ValueError("Shell data must be list of tuples.")

    fList = list()

    for data in dList:

        if not isinstance(data, tuple):

            raise ValueError("Shell dictionary data must be tuple(2).")

        body = data[0]
        labels = data[1]

        lList = [labels] if isinstance(labels, int) else labels

        if not isinstance(lList, list):

            raise ValueError("Face labels must be list.")

        for l in lList:

            if not isinstance(l, int):

                raise ValueError("Face label must be integer.")

            eval("fList.append(body.Shape.Face" + str(l) + ")")

    return Part.Shell(fList)



def makeFaceShell(document, key, shd):

    name = "Shell" + key.capitalize().replace("_", "")
    label = "shell_" + key

    document.recompute()

    shell = document.addObject("Part::Feature", name)
    shell.Label = label
    shell.Shape = faceShell(shd)

    return shell



def makeFuseBody(document, key, bodies):

    name = "Body" + key.capitalize().replace("_", "")
    label = "body_" + key

    # If only one body is given, fuse with itself
    if len(bodies) == 1: bodies.append(bodies[0])

    document.recompute()

    body = document.addObject("Part::MultiFuse", name)
    body.Label = label
    body.Shapes = bodies

    return body



def exportMeshes(obj, dir, prefix, tol=0.1, scale=1.0):

    S = FreeCAD.Matrix()
    S.scale(scale, scale, scale)

    for o in obj:

        mesh = Mesh.Mesh(o.Shape.tessellate(tol))
        #mesh = MeshPart.meshFromShape(Shape=o.Shape, MaxLength=10)

        mesh.transform(S)

        name = prefix + "_" + o.Label + ".stl"

        mesh.write(dir + "/" + name)

# --------------------------------------------------------------------------- #
# --- Classes --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# --- Main ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

def main():

    pass

# --------------------------------------------------------------------------- #

if __name__ == "__main__": main()

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

