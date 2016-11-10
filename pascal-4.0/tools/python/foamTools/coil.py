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
# --- Paths ----------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def pathLoop(data, bdata, bi, p0=None, scale=1.0):
    """

    data-Keys
    ----------
    n : int, Number of bundles
    r : float, Coil radius
    """

    if not 'shape' in bdata:

        raise KeyError("Bundle shape (shape) is missing.")

    if not 'n' in data:

        raise KeyError("Number of edges (n) is missing.")

    if not type(data['n']) ==  int:

        raise KeyError("Number of edges (n) needs to be of type int.")

    if not data['n'] > 0:

        raise ValueError("Number of edges (n) needs to be larger than 0.")

    if not 'r' in data:

        raise KeyError("Coil loop radius (r) is missing.")

    if data['r'] <= 0.0:

        raise ValueError("Coil loop radius (r) must be positive.")

    if not p0: p0 = np.zeros(3)

    rb, zb = bundle[bdata['shape']](bdata, bi)

    points = []

    for i in range(data['n']):

        p = np.zeros(3)

        p[0] = p0[0] + (data['r'] + rb) * m.cos(i*1.0/data['n'] * 2.0*m.pi)
        p[1] = p0[1] + (data['r'] + rb) * m.sin(i*1.0/data['n'] * 2.0*m.pi)
        p[2] = p0[2] + zb

        points.append(scale*p)

    return points

# --------------------------------------------------------------------------- #

path = {'loop': pathLoop}

# --------------------------------------------------------------------------- #
# --- Bundles --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def bundlePointN(data):

    return 1

def bundlePointI(data, I):

    return I

def bundlePoint(data, i, x0=0.0, y0=0.0, scale=1.0):
    """

    data-Keys
    ----------
    """

    return scale*x0, scale*y0

# --------------------------------------------------------------------------- #

def bundleCircleN(data):

    if not 'n' in data:

        raise KeyError("Number of bundles (n) is missing.")

    return data['n']

def bundleCircleI(data, I):

    return I/bundleCircleN(data)

def bundleCircle(data, i, x0=0.0, y0=0.0, scale=1.0):
    """

    data-Keys
    ----------
    n : int, Number of bundles
    r : float, Coil radius
    """

    if not 'n' in data:

        raise KeyError("Number of bundles (n) is missing.")

    if not type(data['n']) ==  int:

        raise KeyError("Number of bundles (n) needs to be of type int.")

    if not data['n'] > 0:

        raise ValueError("Number of bundles (n) needs to be larger than 0.")

    if not 'r' in data:

        raise KeyError("Coil bundle radius (r) is missing.")

    if data['r'] <= 0.0:

        raise ValueError("Coil bundle radius (r) must be positive.")

    if not i < data['n']:

        raise ValueError("Bundle index (i) out of range (max: n).")

    x = x0 + data['r']*m.cos(i*1.0/data['n'] * 2.0*m.pi)
    y = y0 + data['r']*m.sin(i*1.0/data['n'] * 2.0*m.pi)

    return scale*x, scale*y

# --------------------------------------------------------------------------- #

def bundleRectangleN(data, s=4):

    if not 'n' in data:

        raise KeyError("Number of bundles per side (n) is missing.")

    return s*(data['n']-1)

def bundleRectangleI(data, I):

    return I/bundleRectangleN(data)

def bundleRectangle(data, i, x0=0.0, y0=0.0, scale=1.0):
    """

    data-Keys
    ----------
    n : int, Number of bundles per side
    x : float,  Coil size in radial direction
    y : float,  Coil size in axial direction
    """

    if not 'n' in data:

        raise KeyError("Number of bundles per side (n) is missing.")

    if not type(data['n']) ==  int:

        raise KeyError("Number of bundles per side (n) needs to be of type int.")

    if not data['n'] > 1:

        raise ValueError("Number of bundles per side (n) needs to be larger than 1.")

    if not 'x' in data:

        raise KeyError("Coil bundle size (x) is missing.")

    if not 'y' in data:

        raise KeyError("Coil bundle size (y) is missing.")

    if data['x'] <= 0.0 or data['y'] <= 0.0:

        raise ValueError("Coil size (x/y) must be positive.")

    def N(s): return bundleRectangleN(data, s)
    def N0(s): return bundleRectangleN(data, s) + 1

    if not i < N(4):

        raise ValueError("Bundle index (i) out of range (max: 4*(n-1)).")

    x = x0 - data['x']/2.0
    y = y0 - data['y']/2.0

    xi = data['x'] / (data['n']-1)
    yi = data['y'] / (data['n']-1)

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

# --------------------------------------------------------------------------- #

bundleN = {'point': bundlePointN, 'circle': bundleCircleN, 'rectangle': bundleRectangleN}

bundleI = {'point': bundlePointI, 'circle': bundleCircleI, 'rectangle': bundleRectangleI}

bundle = {'point': bundlePoint, 'circle': bundleCircle, 'rectangle': bundleRectangle}

# --------------------------------------------------------------------------- #
# --- Writing --------------------------------------------------------------- #
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

# --------------------------------------------------------------------------- #

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

# --------------------------------------------------------------------------- #

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

