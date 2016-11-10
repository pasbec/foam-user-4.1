#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

import math as m
import numpy as np

from ioInfo import objectIndent, objectHeader, objectFooter

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

def bundleCircleN(n): return n

def bundleCircle(i, n, r, x0=0.0, y0=0.0):

    if not type(n) ==  int:

        raise KeyError("Number of bundles (n) needs to be of type int.")

    if not n > 0:

        raise ValueError("Number of bundles (n) needs to be larger than 0.")

    if r <= 0.0:

        raise ValueError("Coil size (r) must be positive.")

    if not i < n:

        raise ValueError("Bundle index (i) out of range (max: n).")

    x = x0 + r*m.cos(i*1.0/n * 2.0*m.pi)
    y = y0 + r*m.sin(i*1.0/n * 2.0*m.pi)

    return x, y



def bundleRectangleN(n, s=4): return s*(n-1)

def bundleRectangle(i, n, x, y, x0=0.0, y0=0.0):

    if not type(n) ==  int:

        raise KeyError("Number of bundles per side (n) needs to be of type int.")

    if not n > 1:

        raise ValueError("Number of bundles per side (n) needs to be larger than 1.")

    if x <= 0.0 or y <= 0.0:

        raise ValueError("Coil size (x/y) must be positive.")

    def N(s): return bundleRectangleN(n, s)
    def N0(s): return bundleRectangleN(n, s) + 1

    if not i < N(4):

        raise ValueError("Bundle index (i) out of range (max: 4*(n-1)).")

    x = x0 - x/2.0
    y = y0 - y/2.0

    xi = x / (n-1)
    yi = y / (n-1)

    if (i < N0(1)):

        x += i*xi
        y += 0.0

    elif (i > N(1)) and (i < N0(2)):

        x += N(1)*xi
        y += i*yi - N(1)*yi

    elif (i > N(2)) and (i < N0(3)):

        x += N(3)*xi - i*xi
        y += N(1)*yi

    elif (i > N(3)) and (i < N(4)):

        x += 0.0
        y += i*yi - N(3)*yi

    return x, y



bundleN = {'circle': bundleCircleN, 'rectangle': bundleRectangleN}

bundle = {'circle': bundleCircle, 'rectangle': bundleRectangle}

# --------------------------------------------------------------------------- #

def edgeLoopFromPoints(points, off=0):

    edges = list()

    l = len(points) - 1

    for i in range(l):

        io = off + i

        edges.append([io, io+1])

    edges.append([off+l, off])

    return edges



def writeCoilFeatureEdgeMesh(folder, name, points, edges):

    path = folder + '/' + name + '.eMesh'

    with open(path, 'w') as f:

        # Define short indented line with line break
        def ind(iL, cS, eS='\n'): return objectIndent(cS + eS, iLevel=iL)

        f.write(objectHeader(name, 'featureEdgeMesh'))

        f.write(ind(0, '// points:\n'))
        f.write(ind(0, str(len(points))))
        f.write(ind(0, '('))

        for p in points:

            f.write(ind(1, '(' + str(p[0]) + ' ' + str(p[1]) + ' ' + str(p[2])+')'))

        f.write(ind(0, ')\n'))

        f.write(ind(0, '// edges:\n'))
        f.write(ind(0, str(len(edges))))
        f.write(ind(0, '('))

        for e in edges:

            f.write(ind(1, '(' + str(e[0]) + ' '+str(e[1]) + ')'))

        f.write(ind(0, ')\n'))

        f.write(objectFooter())



def writeEdgeBiotSavartProperties(folder, nNonOrth, names, reverse, current, phase):

    name = 'edgeBiotSavartProperties'

    path = folder + '/' + name

    n = len(names)

    if not type(reverse) == list: reverse = [reverse for i in range(n)]
    if not type(current) == list: current = [current for i in range(n)]
    if not type(phase) == list: phase = [phase for i in range(n)]

    with open(path, 'w') as f:

        # Define short indented line with line break
        def ind(iL, cS, eS='\n'): return objectIndent(cS + eS, iLevel=iL)

        f.write(objectHeader(name, 'dictionary'))

        f.write(ind(0, 'nNonOrthogonalCorrectors    ' + str(nNonOrth) + ';\n'))

        f.write(ind(0, 'inductors'))
        f.write(ind(0, '{'))

        for i in range(n):

            def bs(b): return 'true' if b else 'false'

            f.write(ind(1, names[i] + ''))
            f.write(ind(1, '{'))

            f.write(ind(2, 'file       ' + '"' + names[i] + '.eMesh"' + ';'))
            f.write(ind(2, 'reverse    ' + bs(reverse[i]) + ';'))
            f.write(ind(2, 'current    ' + str(current[i]) + ';'))
            f.write(ind(2, 'phase      ' + str(phase[i]) + ';'))

            f.write(ind(1, '}'))

            if not i == n-1: f.write(ind(1, ''))


        f.write(ind(0, '}\n'))

        f.write(objectFooter())



def writeFrequency(folder, value):

    name = 'f0'

    path = folder + '/' + name

    with open(path, 'w') as f:

        # Define short indented line with line break
        def ind(iL, cS, eS='\n'): return objectIndent(cS + eS, iLevel=iL)

        f.write(objectHeader(name, 'uniformDimensionedScalarField'))

        f.write(ind(0, 'dimensions    [0 0 -1 0 0 0 0];\n'))
        f.write(ind(0, 'value         ' + str(value) + ';\n'))

        f.write(objectFooter())

# --------------------------------------------------------------------------- #
# --- End of module --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

