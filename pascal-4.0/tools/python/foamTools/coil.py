#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
# Pascal Beckstein (p.beckstein@hzdr.de)
#
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

def edgeLoopFromPoints(points, off=0):

    edges = list()

    l = len(points) - 1

    for i in range(l):

        io = off + i

        edges.append([io, io+1])

    edges.append([off+l, off])

    return edges

# --------------------------------------------------------------------------- #

def bundlePointN(b):

    return 1

def bundlePointI(b, I):

    return I

def bundlePoint(b, i, x0=0.0, y0=0.0, scale=1.0):

    return scale*x0, scale*y0



def bundleCircleN(b):

    if not 'n' in b:

        raise KeyError("Number of bundles (n) is missing.")

    return b['n']

def bundleCircleI(b, I):

    return I/bundleCircleN(b)

def bundleCircle(b, i, x0=0.0, y0=0.0, scale=1.0):

    if not 'n' in b:

        raise KeyError("Number of bundles (n) is missing.")

    if not type(b['n']) ==  int:

        raise KeyError("Number of bundles (n) needs to be of type int.")

    if not b['n'] > 0:

        raise ValueError("Number of bundles (n) needs to be larger than 0.")

    if not 'r' in b:

        raise KeyError("Coil radius (r) is missing.")

    if b['r'] <= 0.0:

        raise ValueError("Coil radius (r) must be positive.")

    if not i < b['n']:

        raise ValueError("Bundle index (i) out of range (max: n).")

    x = x0 + b['r']*m.cos(i*1.0/b['n'] * 2.0*m.pi)
    y = y0 + b['r']*m.sin(i*1.0/b['n'] * 2.0*m.pi)

    return scale*x, scale*y



def bundleRectangleN(b, s=4):

    if not 'n' in b:

        raise KeyError("Number of bundles per side (n) is missing.")

    return s*(b['n']-1)

def bundleRectangleI(b, I):

    return I/bundleRectangleN(b)

def bundleRectangle(b, i, x0=0.0, y0=0.0, scale=1.0):

    if not 'n' in b:

        raise KeyError("Number of bundles per side (n) is missing.")

    if not type(b['n']) ==  int:

        raise KeyError("Number of bundles per side (n) needs to be of type int.")

    if not b['n'] > 1:

        raise ValueError("Number of bundles per side (n) needs to be larger than 1.")

    if not 'x' in b:

        raise KeyError("Coil size (x) is missing.")

    if not 'y' in b:

        raise KeyError("Coil size (y) is missing.")

    if b['x'] <= 0.0 or b['y'] <= 0.0:

        raise ValueError("Coil size (x/y) must be positive.")

    def N(s): return bundleRectangleN(b, s)
    def N0(s): return bundleRectangleN(b, s) + 1

    if not i < N(4):

        raise ValueError("Bundle index (i) out of range (max: 4*(n-1)).")

    x = x0 - b['x']/2.0
    y = y0 - b['y']/2.0

    xi = b['x'] / (b['n']-1)
    yi = b['y'] / (b['n']-1)

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

    return scale*x, scale*y



bundleN = {'point': bundlePointN, 'circle': bundleCircleN, 'rectangle': bundleRectangleN}

bundleI = {'point': bundlePointI, 'circle': bundleCircleI, 'rectangle': bundleRectangleI}

bundle = {'point': bundlePoint, 'circle': bundleCircle, 'rectangle': bundleRectangle}

# --------------------------------------------------------------------------- #

def writeCoilFeatureEdgeMeshes(case, names, points, edges):

    for n in names.keys():

        writeCoilFeatureEdgeMesh(case, names[n], points[n], edges[n])



def writeCoilFeatureEdgeMesh(case, name, points, edges):

    path = case + '/' + 'constant' + '/' + 'featureEdgeMesh'
    if not os.path.exists(path): os.makedirs(path)
    fullpath = path + '/' + name + '.eMesh'

    with open(fullpath, 'w') as f:

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



def writeEdgeBiotSavartProperties(case, names, b, bs):

    n = len(names)

    if not 'shape' in b:

        raise KeyError("Bundle shape (shape) is missing.")

    if not 'nNonOrth' in bs:

        raise KeyError("Number of non-orthogonal correctors (nNonOrth) is missing.")

    if not 'reverse' in bs:

        raise KeyError("Reverse settings (reverse) are missing.")

    if not 'current' in bs:

        raise KeyError("Current settings (current) are missing.")

    if not 'phase' in bs:

        raise KeyError("Phase settings (phase) are missing.")

    if not len(bs['reverse']) == n:

        raise KeyError("Number of reverse settings is wrong.")

    if not len(bs['current']) == n:

        raise KeyError("Number of current settings is wrong.")

    if not len(bs['phase']) == n:

        raise KeyError("Number of phase settings is wrong.")

    path = case + '/' + 'constant'
    if not os.path.exists(path): os.makedirs(path)
    name = 'edgeBiotSavartProperties'
    fullpath = path + '/' + name

    with open(fullpath, 'w') as f:

        # Define short indented line with line break
        def ind(iL, cS, eS='\n'): return objectIndent(cS + eS, iLevel=iL)

        # Device total current by number of bundles
        def reducedCurrent(I): return bundleI[b['shape']](b, I)

        f.write(objectHeader(name, 'dictionary'))

        f.write(ind(0, 'nNonOrthogonalCorrectors    ' + str(bs['nNonOrth']) + ';\n'))

        f.write(ind(0, 'inductors'))
        f.write(ind(0, '{'))

        for i in range(n):

            def bstr(b): return 'true' if b else 'false'

            f.write(ind(1, names[i] + ''))
            f.write(ind(1, '{'))

            f.write(ind(2, 'file       ' + '"' + names[i] + '.eMesh"' + ';'))
            f.write(ind(2, 'reverse    ' + bstr(bs['reverse'][i]) + ';'))
            f.write(ind(2, 'current    ' + str(reducedCurrent(bs['current'][i])) + ';'))
            f.write(ind(2, 'phase      ' + str(bs['phase'][i]) + ';'))

            f.write(ind(1, '}'))

            if not i == n-1: f.write(ind(1, ''))


        f.write(ind(0, '}\n'))

        f.write(objectFooter())



def writeFrequency(case, value):

    path = case + '/' + 'constant'
    if not os.path.exists(path): os.makedirs(path)
    name = 'f0'
    fullpath = path + '/' + name

    with open(fullpath, 'w') as f:

        # Define short indented line with line break
        def ind(iL, cS, eS='\n'): return objectIndent(cS + eS, iLevel=iL)

        f.write(objectHeader(name, 'uniformDimensionedScalarField'))

        f.write(ind(0, 'dimensions    [0 0 -1 0 0 0 0];\n'))
        f.write(ind(0, 'value         ' + str(value) + ';\n'))

        f.write(objectFooter())

# --------------------------------------------------------------------------- #
# --- End of module --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

