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

def edgeLoopFromPoints(points, start=0):

    edges = list()

    l = len(points) - 1

    for i in range(l):

        I = start + i

        edges.append([I, I + 1])

    edges.append([start + l, start])

    return edges

# --------------------------------------------------------------------------- #
# --- Paths ----------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def pathLoop(pathDict, b, start):
    """

    pathDict-Keys
    ----------
    n : int, Number of edges
    r : float, Coil loop radius
    """

    if not 'n' in pathDict:

        raise KeyError("Number of edges (n) is missing.")

    if not type(pathDict['n']) ==  int:

        raise KeyError("Number of edges (n) needs to be of type int.")

    if not pathDict['n'] > 0:

        raise ValueError("Number of edges (n) needs to be larger than 0.")

    if not 'r' in pathDict:

        raise KeyError("Coil loop radius (r) is missing.")

    if pathDict['r'] <= 0.0:

        raise ValueError("Coil loop radius (r) must be positive.")

    points = list()

    r = b[0] + pathDict['r']
    z = b[1]
    phii = 1.0/pathDict['n'] * 2.0*m.pi

    for i in range(pathDict['n']):

        p = np.zeros(3)

        p[0] = r * m.cos(i * phii)
        p[1] = r * m.sin(i * phii)
        p[2] = z

        points.append(p)

    edges = edgeLoopFromPoints(points, start)

    return points, edges

# --------------------------------------------------------------------------- #

path = {'loop': pathLoop}

# --------------------------------------------------------------------------- #
# --- Bundles --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def bundlePointN(bundleDict):

    return 1

def bundlePointI(bundleDict, I):

    return I

def bundlePoint(bundleDict, i):
    """

    bundleDict-Keys
    ----------
    """

    b = np.zeros(2)

    return b

# --------------------------------------------------------------------------- #

def bundleCircleN(bundleDict):

    if not 'n' in bundleDict:

        raise KeyError("Number of bundles (n) is missing.")

    return bundleDict['n']

def bundleCircleI(bundleDict, I):

    return I/bundleCircleN(bundleDict)

def bundleCircle(bundleDict, i):
    """

    bundleDict-Keys
    ----------
    n : int, Number of bundles
    r : float, Coil bundle radius
    """

    if not 'n' in bundleDict:

        raise KeyError("Number of bundles (n) is missing.")

    if not type(bundleDict['n']) ==  int:

        raise KeyError("Number of bundles (n) needs to be of type int.")

    if not bundleDict['n'] > 0:

        raise ValueError("Number of bundles (n) needs to be larger than 0.")

    if not 'r' in bundleDict:

        raise KeyError("Coil bundle radius (r) is missing.")

    if bundleDict['r'] <= 0.0:

        raise ValueError("Coil bundle radius (r) must be positive.")

    if not i < bundleDict['n']:

        raise ValueError("Coil bundle index (i) out of range (max: n).")

    r = bundleDict['r']
    phii = 1.0/bundleDict['n'] * 2.0*m.pi

    b = r * np.array([m.cos(i * phii), m.sin(i * phii)])

    return b

# --------------------------------------------------------------------------- #

def bundleRectangleN(bundleDict, s=4):

    if not 'n' in bundleDict:

        raise KeyError("Number of bundles per side (n) is missing.")

    return s*(bundleDict['n']-1)

def bundleRectangleI(bundleDict, I):

    return I/bundleRectangleN(bundleDict)

def bundleRectangle(bundleDict, i):
    """

    bundleDict-Keys
    ----------
    n : int, Number of bundles per side
    x : float,  Coil size in radial direction
    y : float,  Coil size in axial direction
    """

    if not 'n' in bundleDict:

        raise KeyError("Number of bundles per side (n) is missing.")

    if not type(bundleDict['n']) ==  int:

        raise KeyError("Number of bundles per side (n) needs to be of type int.")

    if not bundleDict['n'] > 1:

        raise ValueError("Number of bundles per side (n) needs to be larger than 1.")

    if not 'x' in bundleDict:

        raise KeyError("Coil bundle size (x) is missing.")

    if not 'y' in bundleDict:

        raise KeyError("Coil bundle size (y) is missing.")

    if bundleDict['x'] <= 0.0 or bundleDict['y'] <= 0.0:

        raise ValueError("Coil bundle size (x/y) must be positive.")

    def N(s): return bundleRectangleN(bundleDict, s)
    def N0(s): return bundleRectangleN(bundleDict, s) + 1

    if not i < N(4):

        raise ValueError("Coil bundle index (i) out of range (max: 4*(n-1)).")

    b = -0.5 * np.array([bundleDict['x'], bundleDict['y']])

    bi = -2.0 * b / (bundleDict['n'] - 1)

    if (i < N0(1)):

        b[0] += bi[0] * i
        b[1] += bi[1] * 0.0

    elif (i > N(1)) and (i < N0(2)):

        b[0] += bi[0] * N(1)
        b[1] += bi[1] * (i - N(1))

    elif (i > N(2)) and (i < N0(3)):

        b[0] += bi[0] * (N(3) - i)
        b[1] += bi[1] * N(1)

    elif (i > N(3)) and (i < N(4)):

        b[0] += bi[0] * 0.0
        b[1] += bi[1] * (i - N(3))

    return b

# --------------------------------------------------------------------------- #

bundleN = {'point': bundlePointN, 'circle': bundleCircleN, 'rectangle': bundleRectangleN}

bundleI = {'point': bundlePointI, 'circle': bundleCircleI, 'rectangle': bundleRectangleI}

bundle = {'point': bundlePoint, 'circle': bundleCircle, 'rectangle': bundleRectangle}

# --------------------------------------------------------------------------- #
# --- Writing --------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def writeCoilFeatureEdgeMeshes(case, coils):

    for n in coils.keys():

        writeCoilFeatureEdgeMesh(case, coils[n].name,
                                 coils[n].points, coils[n].edges)

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

def writeEdgeBiotSavartProperties(case, coils, nNonOrth=10):

    path = case + '/' + 'constant'
    if not os.path.exists(path): os.makedirs(path)
    name = 'edgeBiotSavartProperties'
    fullpath = path + '/' + name

    with open(fullpath, 'w') as f:

        # Define short indented line with line break
        def ind(iL, cS, eS='\n'): return objectIndent(cS + eS, iLevel=iL)

        # Define write for boolean-strings for OpenFOAM
        def bstr(b): return 'true' if b else 'false'

        f.write(objectHeader(name, 'dictionary'))

        f.write(ind(0, 'nNonOrthogonalCorrectors    ' + str(nNonOrth) + ';\n'))

        f.write(ind(0, 'inductors'))
        f.write(ind(0, '{'))

        for i in range(len(coils)):

            name    = coils[i].name
            reverse = coils[i].reverse
            current = coils[i].bundleCurrent
            phase   = coils[i].phase
            bundles =coils[i].bundles

            f.write(ind(1, name))
            f.write(ind(1, '{'))

            f.write(ind(2, 'file       ' + '"' + name + '.eMesh"' + ';'))
            f.write(ind(2, 'reverse    ' + bstr(reverse) + ';'))
            f.write(ind(2, 'current    ' + str(current) + ';'))
            f.write(ind(2, 'phase      ' + str(phase) + ';'))
            f.write(ind(2, 'bundles    ' + str(bundles) + ';'))

            f.write(ind(1, '}'))

            if not i == len(coils) - 1: f.write(ind(1, ''))


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
# --- Class definitions ----------------------------------------------------- #
# --------------------------------------------------------------------------- #

class coil(object):

    def __init__(self, name, reverse, current, phase,
                 bundleDict, pathDict, **kwargs):

        if not type(name) == str:

            raise KeyError("Name must be of type string")

        if not type(reverse) == bool:
            
            raise KeyError("Reverse switch must be of type bool.")

        if not type(current) in [int, float]:

            raise KeyError("Current must be of type int or float.")

        if not type(phase) in [int, float]:

            raise KeyError("Phase must be of type int or float.")

        self.name = name
        self.reverse = reverse
        self.current = current
        self.bundleCurrent = None
        self.phase = phase
        
        self.bundleDict = bundleDict
        self.pathDict = pathDict

        self.points = list()
        self.edges = list()

        self.compute(**kwargs)

    # ----------------------------------------------------------------------- #

    def compute(self, **kwargs):

        bundleShape = self.bundleDict['shape']
        pathShape   = self.pathDict['shape']

        self.bundles = bundleN[bundleShape](self.bundleDict)
        self.bundleCurrent = bundleI[bundleShape](self.bundleDict, self.current)

        for i in range(self.bundles):

            b = bundle[bundleShape](self.bundleDict, i)

            start = len(self.points)

            p, e = path[pathShape](self.pathDict, b, start)

            self.points += p
            self.edges += e

        self.transform(**kwargs)

    # ----------------------------------------------------------------------- #

    def transform(self, translate=None, rotate=None, scale=None):

        if translate:

            if not type(translate) in [list, np.array]:

                raise KeyError("Translation vector must be a list or an array.")

            if type(translate) is list: translate = np.array(translate)

            if not len(translate) == 3:

                raise KeyError("Translation vector must have exactly 3 components.")

            p = self.points
            for i in range(len(p)): p[i] += translate

        if rotate:

            if not type(rotate) == list and len(rotate) == 2:

                raise KeyError("Rotation data must be a list containing axis and angle.")

            rotAxis  = rotate[0]
            rotAngle = rotate[1]

            if not type(rotAxis) in [list, np.array]:

                raise KeyError("Rotation vector must be a list or an array.")

            if type(rotAxis) is list: rotAxis = np.array(rotAxis)

            if not len(rotAxis) == 3:

                raise KeyError("Rotation vector must have exactly 3 components.")

            from scipy.linalg import expm3, norm

            def R(axis, theta):

                return expm3(np.cross(np.eye(3), theta * axis/norm(axis)))

            rotTheta = rotAngle/180.0 * m.pi
            rotM = R(rotAxis, rotTheta)

            p = self.points
            for i in range(len(p)): p[i] = np.dot(rotM, p[i])

        if scale:

            if not type(scale) in [int, float, list, np.array]:

                raise KeyError("Scalae factor must be of type int, float, list or array.")

            if type(scale) in [int, float]: scale = np.array(3 * [scale])

            if type(scale) is list: scale = np.array(scale)

            if not len(scale) == 3:

                raise KeyError("Scaling vector must have exactly 3 components.")

            p = self.points
            for i in range(len(p)): p[i] *= scale

    # ----------------------------------------------------------------------- #

    def printData(self):

        print("name:", self.name)
        print("reverse:", reverse)
        print("current:", current)
        print("phase:", phase)

        print("pathDict:", self.pathDict)
        print("bundleDict:", self.bundleDict)

        print("points:", self.points)
        print("edges:", self.edges)

# --------------------------------------------------------------------------- #
# --- End of module --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

