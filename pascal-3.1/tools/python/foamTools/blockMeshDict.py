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

        for verticeIndex, verticeLabel in enumerate(self.labels):

            wstr  = "( "
            for component in self.points[verticeIndex]:
                wstr += str(float(component)) + " "
            wstr += ")"
            self.stdout.write(wstr)



class blocks:

    class cBlockTopo:

        origin = [0]

        base = [[0,1],[0,3],[0,4]]

        edgeVertices = [[0,1],[3,2],[7,6],[4,5],
                        [0,3],[1,2],[5,6],[4,7],
                        [0,4],[1,5],[2,6],[3,7]]

        faceVertices = [[0,3,7,4], [1,2,6,5],
                        [0,1,5,4], [2,3,7,6],
                        [0,1,2,3], [4,5,6,7]]

        faceTags = {"x-": 0, "x+": 1,
                    "y-": 2, "y+": 3,
                    "z-": 4, "z+": 5}

    labels = list()
    labelIndex = resizeList()

    topo = cBlockTopo()

    blockVertices = list()
    blockVerticeLabels = list()
    divider = list()
    gradings = list()
    zones = list()

    neighbours = list()
    neighbourFaces = list()
    neighbourEdges = list()

    def __init__(self, stdoutRef, verticesRef):

        self.stdout = stdoutRef
        self.vertices = verticesRef

    def _getSharedVertices(self, vertices, otherVertices):

        shared = list(set(vertices) & set(otherVertices))

        return shared

    def _getNeighbours(self, blockIndex, otherBlockIndex):

        blockVertices = self.blockVertices[blockIndex]
        otherBlockVertices = self.blockVertices[otherBlockIndex]

        neighbourVertices = \
            self._getSharedVertices(blockVertices, otherBlockVertices)

        if len(neighbourVertices) == 4:

            return True, "face", neighbourVertices

        if len(neighbourVertices) == 2:

            return True, "edge", neighbourVertices

        else:

            return False, None, neighbourVertices

    def test(self):

        for blockIndex, blockLabel in enumerate(self.labels):

            print "neighbours: ", blockIndex, blockLabel, self.neighbours[blockIndex]
            print "neighbourFaces: ", blockIndex, blockLabel, self.neighbourFaces[blockIndex]
            print "neighbourEdges: ", blockIndex, blockLabel, self.neighbourEdges[blockIndex]
            print

    def set(self, label, blockVerticeLabels, divider=[1,1,1], grading=[1,1,1], zone=""):

        # Convert blockVerticeLabels in blockVertices
        blockVertices = [ self.vertices.labelIndex[l] for l in blockVerticeLabels ]

        try:

            blockIndex = self.labels.index(label)

            self.blockVertices[blockIndex] = blockVertices
            self.blockVerticeLabels[blockIndex] = blockVerticeLabels
            self.divider[blockIndex] = divider
            self.gradings[blockIndex] = grading
            self.zones[blockIndex] = zone
            self.neighbours[blockIndex] = None
            self.neighbourFaces[blockIndex] = [None for i in range(6)]
            self.neighbourEdges[blockIndex] = [None for i in range(12)]

        except:

            blockIndex = len(self.labels)
            self.labels.append(label)
            self.labelIndex[label] = blockIndex

            self.blockVertices.append(blockVertices)
            self.blockVerticeLabels.append(blockVerticeLabels)
            self.divider.append(divider)
            self.gradings.append(grading)
            self.zones.append(zone)
            self.neighbours.append(None)
            self.neighbourFaces.append([None for i in range(6)])
            self.neighbourEdges.append([None for i in range(12)])

        # Identify all neighbours of current block
        for otherBlockIndex, otherBlockLabel in enumerate(self.labels):

            if not blockIndex == otherBlockIndex:

                neighbours, neighbourType, neighbourVertices = \
                    self._getNeighbours(blockIndex, otherBlockIndex)

                if neighbourType == "face":

                    try:
                        self.neighbours[blockIndex].index(otherBlockIndex)
                    except:
                        if self.neighbours[blockIndex] == None:
                            self.neighbours[blockIndex] = []
                        self.neighbours[blockIndex].append(otherBlockIndex)
                    try:
                        self.neighbours[otherBlockIndex].index(blockIndex)
                    except:
                        if self.neighbours[otherBlockIndex] == None:
                            self.neighbours[otherBlockIndex] = []
                        self.neighbours[otherBlockIndex].append(blockIndex)

                    # Identify all face neighbours
                    blockVertices = self.blockVertices[blockIndex]
                    otherBlockVertices = self.blockVertices[otherBlockIndex]

                    for localFaceIndex, localFace in enumerate(self.topo.faceVertices):

                        globalFace = \
                            [ blockVertices[i] for i in localFace ]
                        otherGlobalFace = \
                            [ otherBlockVertices[i] for i in localFace ]

                        sharedFaceVertices = \
                            self._getSharedVertices(neighbourVertices, globalFace)
                        otherSharedFaceVertices = \
                            self._getSharedVertices(neighbourVertices, otherGlobalFace)

                        if len(sharedFaceVertices) == 4:
                            self.neighbourFaces[blockIndex][localFaceIndex] = otherBlockIndex
                        if len(otherSharedFaceVertices) == 4:
                            self.neighbourFaces[otherBlockIndex][localFaceIndex] = blockIndex

                    # Identify all edge neighbours
                    blockVertices = self.blockVertices[blockIndex]
                    otherBlockVertices = self.blockVertices[otherBlockIndex]

                    for localEdgeIndex, localEdge in enumerate(self.topo.edgeVertices):

                        globalEdge = \
                            [ blockVertices[i] for i in localEdge ]
                        otherGlobalEdge = \
                            [ otherBlockVertices[i] for i in localEdge ]
                        
                        sharedEdgeVertices = \
                            self._getSharedVertices(neighbourVertices, globalEdge)
                        otherSharedEdgeVertices = \
                            self._getSharedVertices(neighbourVertices, otherGlobalEdge)

                        if len(sharedEdgeVertices) == 2:
                            try:
                                self.neighbourEdges[blockIndex][localEdgeIndex].index(otherBlockIndex)
                            except:
                                if self.neighbourEdges[blockIndex][localEdgeIndex] == None:
                                    self.neighbourEdges[blockIndex][localEdgeIndex] = []
                                self.neighbourEdges[blockIndex][localEdgeIndex].append(otherBlockIndex)
                        if len(otherSharedEdgeVertices) == 2:
                            try:
                                self.neighbourEdges[otherBlockIndex][localEdgeIndex].index(blockIndex)
                            except:
                                if self.neighbourEdges[otherBlockIndex][localEdgeIndex] == None:
                                    self.neighbourEdges[otherBlockIndex][localEdgeIndex] = []
                                self.neighbourEdges[otherBlockIndex][localEdgeIndex].append(blockIndex)

# TODO [High]: References???

    def write(self):

        def write(string): sys.stdout.write(string)

        for blockIndex, blockLabel in enumerate(self.labels):

            self.stdout.write("hex", end="")

            write(" ")

            wstr  = "( "
            for vertice in self.blockVertices[blockIndex]:
                wstr += str(vertice) + " "
            wstr += ")"
            write(wstr)

            write(" ")

            zone = self.zones[blockIndex]
            if zone:
                write(str(zone) + " ")

            wstr  = "( "
            for divide in self.divider[blockIndex]:
                wstr += str(divide) + " "
            wstr += ")"
            write(wstr)

            write(" ")

            if len(self.gradings[blockIndex]) == 3:
                write("simpleGrading")
            elif len(self.gradings[blockIndex]) == 12:
                write("edgeGrading")

            write(" ")

            wstr  = "( "
            for grading in self.gradings[blockIndex]:
                wstr += str(grading) + " "
            wstr += ")"
            self.stdout.write(wstr, ind=False)



class boundaryFaces:


    labels = list()
    labelIndex = resizeList()

    faces = list()
    faceLabels = list()
    boundary = list()

    def __init__(self, stdoutRef, verticesRef, blocksRef):

        self.stdout = stdoutRef
        self.vertices = verticesRef
        self.blocks = blocksRef

    def set(self, boundary, par1, par2=None):

        if type(par1) == int and par2:

            blockLabel = par1
            blockFaceTag = par2

            blockFaceIndex = self.blocks.topo.faceTags[blockFaceTag]

            blockLocalFace = self.blocks.topo.faceVertices[blockFaceIndex]
            blockIndex = self.blocks.labelIndex[blockLabel]
            blockVertices = self.blocks.blockVertices[blockIndex]
            blockGlobalFace = [ blockVertices[i] for i in blockLocalFace ]

            blockGlobalFaceList = [ blockGlobalFace ]

        elif type(par1) == list and par2:

            blockLabelList = par1
            blockFaceTag = par2

            blockFaceIndex = self.blocks.topo.faceTags[blockFaceTag]

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

            # Convert blockGlobalFace in blockFaceLabels
            blockFaceLabels = [ self.vertices.labelIndex[l] for l in blockGlobalFace ]

            #try: # TODO Label was removed from parameters

                #faceIndex = self.labels.index(label)

                #self.boundary[faceIndex] = boundary
                #self.faces[faceIndex] = blockGlobalFace

            #except:

            label = len(self.labels)

            faceIndex = len(self.labels)
            self.labels.append(label)
            self.labelIndex[label] = faceIndex

            self.boundary.append(boundary)
            self.faces.append(blockGlobalFace)
            self.faceLabels.append(blockFaceLabels)

    def write(self):

        self.stdout.write("faces")
        self.stdout.write("(")
        self.stdout.indentLevel += 1

        for faceIndex, faceLabel in enumerate(self.labels):

            if self.boundary[faceIndex] == self.stdout.subdict.boundary:

                wstr  = "( "
                for vertice in self.faces[faceIndex]:
                    wstr += str(vertice) + " "
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

