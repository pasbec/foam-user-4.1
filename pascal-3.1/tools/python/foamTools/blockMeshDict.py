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

class cResizeList(list):

  def __init__(self, default):

    self._default = default

  def __setitem__(self, key, value):

    if key >= len(self):
      self += [self._default] * (key - len(self) + 1)
    super(cResizeList, self).__setitem__(key, value)



class cStdoutBase:

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



class cVertices:

    labels = []
    labelIndex  = cResizeList(None)

    points = []

    def __init__(self, stdoutRef):

        self.stdout = stdoutRef

    def set(self, label, point):
    
        if label >= len(self.labels):
    
            self.labelIndex[label] = len(self.labels)
            self.labels.append(label)
    
            self.points.append(point)
    
        else:
    
            index = labelIndex[label]
            points[index] = point

    def write(self):
    
        for index, label in enumerate(self.labels):
    
            wstr  = "( "
            for component in self.points[index]:
                wstr += str(float(component)) + " "
            wstr += ")"
            self.stdout.write(wstr)



class cBlocks:

    class cBlockInfo:
    
        faceVertices = {0: [0,3,7,4],
                        1: [1,2,6,5],
                        2: [0,1,5,4],
                        3: [2,3,7,6],
                        4: [0,1,2,3],
                        5: [4,5,6,7]}
    
    
        faceTags = {"x-": 0,
                    "x+": 1,
                    "y-": 2,
                    "y+": 3,
                    "z-": 4,
                    "z+": 5}


    labels = []
    labelIndex = cResizeList(None)

    blockVertices = []
    divider = []
    gradings = []
    zones = []

    def __init__(self, stdoutRef, verticesRef):

        self.info = self.cBlockInfo()
        self.stdout = stdoutRef
        self.vertices = verticesRef

    def set(self, label, blockVertices, divider, grading=[1,1,1], zone=""):
    
        if label >= len(self.labels):
    
            self.labelIndex[label] = len(self.labels)
            self.labels.append(label)

            self.blockVertices.append(blockVertices)
            self.divider.append(divider)
            self.gradings.append(grading)
            self.zones.append(zone)
    
        else:

            index = labelIndex[label]
    
            blockVertices[index] = blockVertices
            divider[index] = divider
            gradings[index] = grading
            zones[index] = zone

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



class cBlockMeshDict:

    class cStdout(cStdoutBase):

        class cSubDict:

            opened = False
            boundary = False

        def __init__(self):

            self.subdict = self.cSubDict()

        def boundarySubDictHeader(self, name, typename="empty"):

            self.subdict.boundary = True
            self.write(str(name))
            self.write("{")
            self.indentLevel += 1
            self.write("type" + " " + str(typename) + ";")
            self.write("faces")
            self.write("(")
            self.indentLevel += 1

        def boundarySubDictFooter(self, line = True):

            self.indentLevel -= 1
            self.write(");")
            self.indentLevel -= 1
            self.write("}")
            self.subdict.boundary = False
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

        self.stdout = self.cStdout()
        self.vertices = cVertices(self.stdout)
        self.blocks = cBlocks(self.stdout, self.vertices)

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
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

def printFace(d):

    if len(d) == 2:

        block = d[0]
        plane = d[1]
        
    if len(d) == 4:

        faceVertices = d

        s  = "( "
        for vertice in faceVertices:
            s += str(verticeMap[vertice]) + " "
        s += ")"
        printi(s)



# --------------------------------------------------------------------------- #
# --- End of module --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

