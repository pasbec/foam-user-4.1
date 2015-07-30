#!/usr/bin/python
# Module with helper functions to write a blockdict with python
# July 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

from foamTools.ioInfo import objectIndent, objectHeader, objectFooter

# --------------------------------------------------------------------------- #
# --- Class definitions ----------------------------------------------------- #
# --------------------------------------------------------------------------- #

class resizeList(list):

  def __init__(self, default=None):

    self.default = default

  def __setitem__(self, key, value):

    if key >= len(self):
      self += [self.default] * (key - len(self) + 1)
    super(resizeList, self).__setitem__(key, value)



class stdoutBase:

    indentLevel = 0

    def line(self):
        sys.stdout.write("\n")

    def indent(self, level, string, end):

        return objectIndent(string + end, iLevel=level)

    def write(self, string, ind=True, end="\n"):

        level = self.indentLevel
        if not ind: level = 0

        wstr = self.indent(level,str(string), end)
        sys.stdout.write(wstr)



class vertices:

    labels = list()
    labelIndex  = resizeList()

    points = list()

    def __init__(self, stdoutRef):

        self.stdout = stdoutRef

    def set(self, label, point):

        try:

            verticeIndex = self.labels.index(label)

            self.points[verticeIndex] = point

        except:

            verticeIndex = len(self.labels)
            self.labels.append(label)
            self.labelIndex[label] = verticeIndex

            self.points.append(point)

    def write(self):

        for index, label in enumerate(self.labels):

            wstr  = "( "
            for component in self.points[index]:
                wstr += str(float(component)) + " "
            wstr += ")"
            self.stdout.write(wstr)



class blocks:

    class cBlockTopo:

        faceVertices = {"x-": [0,3,7,4],
                        "x+": [1,2,6,5],
                        "y-": [0,1,5,4],
                        "y+": [2,3,7,6],
                        "z-": [0,1,2,3],
                        "z+": [4,5,6,7]}

    labels = list()
    labelIndex = resizeList()

    topo = cBlockTopo()

    blockVertices = list()
    divider = list()
    gradings = list()
    zones = list()

    def __init__(self, stdoutRef, verticesRef):

        self.stdout = stdoutRef
        self.vertices = verticesRef

    def set(self, label, blockVertices, divider=[1,1,1], grading=[1,1,1], zone=""):

        try:

            blockIndex = self.labels.index(label)

            self.blockVertices[blockIndex] = blockVertices
            self.divider[blockIndex] = divider
            self.gradings[blockIndex] = grading
            self.zones[blockIndex] = zone

        except:

            blockIndex = len(self.labels)
            self.labels.append(label)
            self.labelIndex[label] = blockIndex

            self.blockVertices.append(blockVertices)
            self.divider.append(divider)
            self.gradings.append(grading)
            self.zones.append(zone)

    def write(self):

        def write(string): sys.stdout.write(string)

        for index, label in enumerate(self.labels):

            self.stdout.write("hex", end="")

            write(" ")

            wstr  = "( "
            for vertice in self.blockVertices[index]:
                wstr += str(self.vertices.labelIndex[vertice]) + " "
            wstr += ")"
            write(wstr)

            write(" ")

            zone = self.zones[index]
            if zone:
                write(str(zone) + " ")

            wstr  = "( "
            for divide in self.divider[index]:
                wstr += str(divide) + " "
            wstr += ")"
            write(wstr)

            write(" ")

            if len(self.gradings[index]) == 3:
                write("simpleGrading")
            elif len(self.gradings[index]) == 12:
                write("edgeGrading")

            write(" ")

            wstr  = "( "
            for grading in self.gradings[index]:
                wstr += str(grading) + " "
            wstr += ")"
            self.stdout.write(wstr, ind=False)



class boundaryFaces:


    labels = list()
    labelIndex = resizeList()

    faces = list()
    boundary = list()

    def __init__(self, stdoutRef, verticesRef, blocksRef):

        self.stdout = stdoutRef
        self.vertices = verticesRef
        self.blocks = blocksRef

    def set(self, boundary, par1, par2=None):

        if type(par1) == int and par2:

            blockLabel = par1
            blockFaceIndex = par2

            blockLocalFace = self.blocks.topo.faceVertices[blockFaceIndex]
            blockIndex = self.blocks.labelIndex[blockLabel]
            blockVertices = self.blocks.blockVertices[blockIndex]
            blockGlobalFace = [ blockVertices[i] for i in blockLocalFace ]

            blockGlobalFaceList = [ blockGlobalFace ]

        elif type(par1) == list and par2:

            blockLabelList = par1
            blockFaceIndex = par2

            blockGlobalFaceList = []

            for blockLabel in blockLabelList:

                blockLocalFace = self.blocks.topo.faceVertices[blockFaceIndex]
                blockIndex = self.blocks.labelIndex[blockLabel]
                blockVertices = self.blocks.blockVertices[blockIndex]
                blockGlobalFace = [ blockVertices[i] for i in blockLocalFace ]

                blockGlobalFaceList.append(blockGlobalFace)

        elif type(par1) == list and len(par1) == 4 and not par2:

            blockGlobalFace = par1

            blockGlobalFaceList = [ blockGlobalFace ]

        else:

            raise KeyError()

        for blockGlobalFace in blockGlobalFaceList:

            #try: # TODO Label was removed from parameters

                #faceIndex = self.labels.index(label)

                #self.boundary[faceIndex] = boundary
                #self.faces[faceIndex] = blockGlobalFace

            #except:

            label = len(self.labels) # TODO

            faceIndex = len(self.labels)
            self.labels.append(label)
            self.labelIndex[label] = faceIndex

            self.boundary.append(boundary)
            self.faces.append(blockGlobalFace)

    def write(self):

        self.stdout.write("faces")
        self.stdout.write("(")
        self.stdout.indentLevel += 1

        for index, label in enumerate(self.labels):

            if self.boundary[index] == self.stdout.subdict.boundary:

                wstr  = "( "
                for vertice in self.faces[index]:
                    wstr += str(self.vertices.labelIndex[vertice]) + " "
                wstr += ")"
                self.stdout.write(wstr)

        self.stdout.indentLevel -= 1
        self.stdout.write(");")



class blockMeshDict:

    class stdoutDict(stdoutBase):

        class subDictData:

            opened = False
            boundary = None

        subdict = subDictData()

        def boundarySubDictHeader(self, name, typename="empty"):

            self.subdict.boundary = name
            self.write(str(name))
            self.write("{")
            self.indentLevel += 1
            self.write("type" + " " + str(typename) + ";")

        def boundarySubDictFooter(self, line = True):

            self.indentLevel -= 1
            self.write("}")
            self.subdict.boundary = None
            if line: self.line()

        def boundarySubDict(self, name, typename="empty"):

            if self.subdict.boundary: self.boundarySubDictFooter()
            self.boundarySubDictHeader(name, typename)
            return True

        def subDictHeader(self, name):

            self.subdict.opened = True
            self.write(str(name))
            self.write("(")
            self.indentLevel += 1

        def subDictFooter(self, line = True):

            self.indentLevel -= 1
            self.write(");")
            self.subdict.opened = False
            if line: self.line()

        def subDict(self, name):

            if self.subdict.boundary: self.boundarySubDictFooter(line=False)
            if self.subdict.opened: self.subDictFooter()
            self.subDictHeader(name)
            return True

        def header(self, scale):

            wstr = objectHeader("blockMeshDict", "dictionary")
            sys.stdout.write(wstr)
            wstr = "convertToMeters" + " " + str(scale) + ";"
            self.write(wstr)
            self.line()

        def footer(self):

            if self.subdict.opened: self.subDictFooter()
            wstr = objectFooter()
            sys.stdout.write(wstr)

    def __init__(self):

        self.stdout = self.stdoutDict()
        self.vertices = vertices(self.stdout)
        self.blocks = blocks(self.stdout, self.vertices)
        self.boundaryFaces = boundaryFaces(self.stdout, self.vertices, self.blocks)

    def boundarySubDict(self, name, typename):

        self.stdout.boundarySubDict(name, typename)
        return True

    def subDict(self, name):

        self.stdout.subDict(name)
        return True

    def header(self, scale):

        self.stdout.header(scale)

    def footer(self):

        self.stdout.footer()

# --------------------------------------------------------------------------- #
# --- End of module --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

