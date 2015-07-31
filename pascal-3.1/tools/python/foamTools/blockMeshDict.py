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

    # ----------------------------------------------------------------------- #

    def __setitem__(self, key, value):

        if key >= len(self):

            self += [self.default] * (key - len(self) + 1)

        super(resizeList, self).__setitem__(key, value)



# --------------------------------------------------------------------------- #



class stdoutBase:

    indentLevel = 0

    def line(self):
        sys.stdout.write("\n")

    # ----------------------------------------------------------------------- #

    def indent(self, level, string, end):

        if not (type(level) == int \
            and type(string) == str \
            and type(end) == str):

            raise KeyError()

        return objectIndent(string + end, iLevel=level)

    # ----------------------------------------------------------------------- #

    def write(self, string, ind=True, end="\n"):

        if not (type(string) == str \
            and type(ind) == bool \
            and type(end) == str):

            raise KeyError()

        level = self.indentLevel
        if not ind: level = 0

        wstr = self.indent(level,str(string), end)
        sys.stdout.write(wstr)



# --------------------------------------------------------------------------- #



class vertices:

    labels = list()
    labelIndex  = resizeList()

    points = list()

    # ----------------------------------------------------------------------- #

    def __init__(self, stdoutRef):

        self.stdout = stdoutRef

    # ----------------------------------------------------------------------- #

    def set(self, label, point):

        if not (type(label) == int \
            and type(point) == list \
            and len(point) == 3):

            raise KeyError()

        # Vertice already exists
        try:

            verticeIndex = self.labels.index(label)

            self.points[verticeIndex] = point

        # New vertice
        except:

            verticeIndex = len(self.labels)
            self.labels.append(label)
            self.labelIndex[label] = verticeIndex

            self.points.append(point)

    # ----------------------------------------------------------------------- #

    def write(self):

        for verticeIndex, verticeLabel in enumerate(self.labels):

            wstr  = "( "
            for component in self.points[verticeIndex]:
                wstr += str(float(component)) + " "
            wstr += ")"
            self.stdout.write(wstr)



# --------------------------------------------------------------------------- #



class blocks:

    class cBlockTopo:

        originVertice = 0

        base = [i for i in range(3)]

        baseVertices = [[0,1],[0,3],[0,4]]

        baseTags = {"x": 0, "y": 1, "z": 2}

        baseGetTags = {0: "x", 1: "y", 2: "z"}

        edges = [i for i in range(12)]

        edgeVertices = [[0,1],[3,2],[7,6],[4,5],
                        [0,3],[1,2],[5,6],[4,7],
                        [0,4],[1,5],[2,6],[3,7]]

        faces = [i for i in range(6)]

        faceVertices = [[0,3,7,4], [1,2,6,5],
                        [0,1,5,4], [2,3,7,6],
                        [0,1,2,3], [4,5,6,7]]

        faceTags = {"x-": 0, "x+": 1,
                    "y-": 2, "y+": 3,
                    "z-": 4, "z+": 5}

        faceGetTags = {0: "x-", 1: "x+",
                       2: "y-", 3: "y+",
                       4: "z-", 5: "z+"}

        faceBaseTags = {"x": [0,1], "y": [2,3], "z": [4,5]}

        faceGetBaseTags = {0: "x", 1: "x",
                           2: "y", 3: "y",
                           4: "z", 5: "z"}

        faceOpposites = {0: 1, 1: 0,
                         2: 3, 3: 2,
                         4: 5, 5: 4}

        faceBase = {0: 0, 1: 0,
                    2: 1, 3: 1,
                    4: 2, 5: 2}

    # ----------------------------------------------------------------------- #

    labels = list()
    labelIndex = resizeList()

    topo = cBlockTopo()

    blockVertices = list()
    blockVerticeLabels = list()
    divider = list()
    gradings = list()
    zones = list()

    neighbours = list()
    faceNeighbours = list()
    edgeNeighbours = list()

    # ----------------------------------------------------------------------- #

    def __init__(self, stdoutRef, verticesRef):

        self.stdout = stdoutRef
        self.vertices = verticesRef

    # ----------------------------------------------------------------------- #

    def _getSharedVertices(self, vertices, otherVertices):

        if not (type(vertices) == list \
            and type(otherVertices) == list):

            raise KeyError()

        return list(set(vertices) & set(otherVertices))

    # ----------------------------------------------------------------------- #

    def _getNeighbours(self, blockIndex, otherBlockIndex):

        if not (type(blockIndex) == int \
            and type(otherBlockIndex) == int):

            raise KeyError()

        # Get both sets of block vertices
        blockVertices = self.blockVertices[blockIndex]
        otherBlockVertices = self.blockVertices[otherBlockIndex]

        # Find shared vertices
        neighbourVertices = \
            self._getSharedVertices(blockVertices, otherBlockVertices)

        # Face neighbours
        if len(neighbourVertices) == 4:

            return True, "face", neighbourVertices

        # Edge neighbours
        if len(neighbourVertices) == 2:

            return True, "edge", neighbourVertices

        # No neighbours
        else:

            return False, None, neighbourVertices

    # ----------------------------------------------------------------------- #

# TODO [High]: Split "_setNeighbourData" into smaller pieces

    def _setNeighbourData(self, blockIndex):

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
                            self.faceNeighbours[blockIndex][localFaceIndex] = otherBlockIndex
                        if len(otherSharedFaceVertices) == 4:
                            self.faceNeighbours[otherBlockIndex][localFaceIndex] = blockIndex

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
                                self.edgeNeighbours[blockIndex][localEdgeIndex].index(otherBlockIndex)
                            except:
                                if self.edgeNeighbours[blockIndex][localEdgeIndex] == None:
                                    self.edgeNeighbours[blockIndex][localEdgeIndex] = []
                                self.edgeNeighbours[blockIndex][localEdgeIndex].append(otherBlockIndex)
                        if len(otherSharedEdgeVertices) == 2:
                            try:
                                self.edgeNeighbours[otherBlockIndex][localEdgeIndex].index(blockIndex)
                            except:
                                if self.edgeNeighbours[otherBlockIndex][localEdgeIndex] == None:
                                    self.edgeNeighbours[otherBlockIndex][localEdgeIndex] = []
                                self.edgeNeighbours[otherBlockIndex][localEdgeIndex].append(blockIndex)

    # ----------------------------------------------------------------------- #

    def set(self, label, verticeLabels, divider=None, grading=None, zone=None):

        if not (type(label) == int \
            and type(verticeLabels) == list \
            and len(verticeLabels) == 8
            and (divider == None or \
                (type(divider) == list and len(divider) == 3)) \
            and (grading == None or \
                (type(grading) == list \
                and (len(grading) == 3 or len(grading) == 12))) \
            and (zone == None or type(zone) == str)):

            raise KeyError()

        blockVertices = [ self.vertices.labelIndex[l] for l in verticeLabels ]

        if not divider: divider = [1,1,1]

        if not grading: grading = [1.0,1.0,1.0]

        # ------------------------------------------------------------------- #

        # Block already exists
        try:

            blockIndex = self.labels.index(label)

            self.blockVertices[blockIndex] = blockVertices
            self.blockVerticeLabels[blockIndex] = verticeLabels
            self.divider[blockIndex] = divider
            self.gradings[blockIndex] = grading
            self.zones[blockIndex] = zone
            self.neighbours[blockIndex] = None
            self.faceNeighbours[blockIndex] = [None for i in range(6)]
            self.edgeNeighbours[blockIndex] = [None for i in range(12)]

        # New block
        except:

            blockIndex = len(self.labels)
            self.labels.append(label)
            self.labelIndex[label] = blockIndex

            self.blockVertices.append(blockVertices)
            self.blockVerticeLabels.append(verticeLabels)
            self.divider.append(divider)
            self.gradings.append(grading)
            self.zones.append(zone)
            self.neighbours.append(None)
            self.faceNeighbours.append([None for i in range(6)])
            self.edgeNeighbours.append([None for i in range(12)])

        self._setNeighbourData(blockIndex)

    # ----------------------------------------------------------------------- #

# TODO [High]: Split "setDivider" into smaller pieces

    def setDivider(self, blockLabels, par1, par2=None):

        if type(blockLabels) == int: blockLabels = [ blockLabels ]

        if not type(blockLabels) == list: raise KeyError()

        blockIndices = [ self.labelIndex[l] for l in blockLabels ]

        for blockIndex in blockIndices:

            # Divider given
            if not par2:

                if type(par1) == list:

                    blockDivider = par1

                    # Set divider for current block
                    self.divider[blockIndex] = blockDivider

                else: raise KeyError()

            # Direction and divide for this direction
            else:

                # Direction as base
                if type(par1) == int:

                    blockDivideBase = par1

                # Direction as baseTag
                elif type(par1) == str:

                    blockDivideBaseTag = par1

                    blockDivideBase = self.topo.baseTags[blockDivideBaseTag]

                else: raise KeyError()

                # Divide
                if type(par2) == int:

                    blockDivide = par2

                else: raise KeyError()

                # Set divide for current block in given direction
                self.divider[blockIndex][blockDivideBase] = blockDivide

                # Get block divider from current block
                blockDivider = self.divider[blockIndex]

        # ------------------------------------------------------------------- #

            # Seed blocks
            seedBlockIndices = []

            # Sync divider in all directions
            for faceIndex in self.topo.faces:

                # Get current and perpendicular directions
                faceBase = self.topo.faceBase[faceIndex]
                faceBasePerpendicular = \
                    [ i for i in self.topo.base if not i == faceBase ]
                faceBaseP0 = faceBasePerpendicular[0]
                faceBaseP1 = faceBasePerpendicular[1]

                # Extract divider in othogonal directions
                blockDivideP0 = blockDivider[faceBaseP0]
                blockDivideP1 = blockDivider[faceBaseP1]

                thisFaceIndex = faceIndex
                thisBlockIndex = blockIndex

                # Hop in current direction given by faceBase
                while True:

                    # Set divide for current block in perpendicular directions
                    # of current direction
                    self.divider[thisBlockIndex][faceBaseP0] = blockDivideP0
                    self.divider[thisBlockIndex][faceBaseP1] = blockDivideP1

                    # Get next neighbour and exit if there is none
                    nextBlockIndex = \
                        self.faceNeighbours[thisBlockIndex][thisFaceIndex]

                    if nextBlockIndex == None: break

                    # Record seed blocks
                    seedBlockIndices.append(nextBlockIndex)

                    # Update indices for next block to become current
                    thisFaceIndex = \
                        self.faceNeighbours[nextBlockIndex].index(thisBlockIndex)
                    thisFaceIndex = \
                        self.topo.faceOpposites[thisFaceIndex]
                    thisBlockIndex = nextBlockIndex

             # Sync divider in all planes (dircetions for all seed blocks)
            for seedBlockIndex in seedBlockIndices:

                seedBlockDivider = self.divider[seedBlockIndex]

                # Sync divider in all directions
                for faceIndex in self.topo.faces:

                    # Get current and perpendicular directions
                    faceBase = self.topo.faceBase[faceIndex]
                    faceBasePerpendicular = \
                        [ i for i in self.topo.base if not i == faceBase ]
                    faceBaseP0 = faceBasePerpendicular[0]
                    faceBaseP1 = faceBasePerpendicular[1]

                    # Extract divider in othogonal directions
                    seedBlockDivideP0 = seedBlockDivider[faceBaseP0]
                    seedBlockDivideP1 = seedBlockDivider[faceBaseP1]

                    thisFaceIndex = faceIndex
                    thisBlockIndex = seedBlockIndex

                    # Hop in current direction given by faceBase
                    while True:

                        # Set divide for current block in perpendicular directions
                        # of current direction
                        self.divider[thisBlockIndex][faceBaseP0] = seedBlockDivideP0
                        self.divider[thisBlockIndex][faceBaseP1] = seedBlockDivideP1

                        # Get next neighbour and exit if there is none
                        nextBlockIndex = \
                            self.faceNeighbours[thisBlockIndex][thisFaceIndex]

                        if nextBlockIndex == None: break

                        # Update indices for next block to become current
                        thisFaceIndex = \
                            self.faceNeighbours[nextBlockIndex].index(thisBlockIndex)
                        thisFaceIndex = \
                            self.topo.faceOpposites[thisFaceIndex]
                        thisBlockIndex = nextBlockIndex

    # ----------------------------------------------------------------------- #

    def write(self):

        def write(string): sys.stdout.write(string)

        for blockIndex, blockLabel in enumerate(self.labels):

            self.stdout.write("hex", end="")

            write(" ")

        # ------------------------------------------------------------------- #

            wstr  = "( "
            for vertice in self.blockVertices[blockIndex]:
                wstr += str(vertice) + " "
            wstr += ")"
            write(wstr)

            write(" ")

        # ------------------------------------------------------------------- #

            zone = self.zones[blockIndex]
            if zone:
                write(str(zone) + " ")

        # ------------------------------------------------------------------- #

            wstr  = "( "
            for divide in self.divider[blockIndex]:
                wstr += str(divide) + " "
            wstr += ")"
            write(wstr)

            write(" ")

        # ------------------------------------------------------------------- #

            if len(self.gradings[blockIndex]) == 3:
                write("simpleGrading")
            elif len(self.gradings[blockIndex]) == 12:
                write("edgeGrading")

            write(" ")

        # ------------------------------------------------------------------- #

            wstr  = "( "
            for grading in self.gradings[blockIndex]:
                wstr += str(grading) + " "
            wstr += ")"
            self.stdout.write(wstr, ind=False)

    # ----------------------------------------------------------------------- #

    def printData(self, blockLabels=None):

        # Only print data for given block
        if blockLabels:

            if type(blockLabels) == int: blockLabels = [ blockLabels ]

            if not type(blockLabels) == list: raise KeyError()

            blockIndices = [ self.labelIndex[l] for l in blockLabels ]

            localBlockIndices = blockIndices

        # Print all blocks
        else:

            localBlockIndices = [ self.labelIndex[l] for l in self.labels ]

        # ------------------------------------------------------------------- #

        for blockIndex in localBlockIndices:

            print "index:", blockIndex
            print "label:", self.labels[blockIndex]
            print "vertices: ", self.blockVertices[blockIndex]
            print "verticeLabels: ", self.blockVerticeLabels[blockIndex]
            print "divider: ", self.divider[blockIndex]
            print "gradings: ", self.gradings[blockIndex]
            print "zone: ", self.zones[blockIndex]
            print "neighbours: ", self.neighbours[blockIndex]
            print "faceNeighbours: ", self.faceNeighbours[blockIndex]
            print "edgeNeighbours: ", self.edgeNeighbours[blockIndex]
            print



# --------------------------------------------------------------------------- #



class boundaryFaces:


    labels = list()
    labelIndex = resizeList()

    faces = list()
    faceLabels = list()
    boundary = list()

    # ----------------------------------------------------------------------- #

    def __init__(self, stdoutRef, verticesRef, blocksRef):

        self.stdout = stdoutRef
        self.vertices = verticesRef
        self.blocks = blocksRef

    # ----------------------------------------------------------------------- #

    def set(self, boundary, par1, par2=None):

        if not type(boundary) == str: raise KeyError()

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

        else: raise KeyError()

        # ------------------------------------------------------------------- #

        for blockGlobalFace in blockGlobalFaceList:

            # Convert blockGlobalFace in blockFaceLabels
            blockFaceLabels = [ self.vertices.labelIndex[l] for l in blockGlobalFace ]

# TODO [Low]: Reorganize! Label has been removed from parameters (not needed?)
            #try:

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

    # ----------------------------------------------------------------------- #

    def write(self):

        self.stdout.write("faces")
        self.stdout.write("(")
        self.stdout.indentLevel += 1

        # ------------------------------------------------------------------- #

        for faceIndex, faceLabel in enumerate(self.labels):

            if self.boundary[faceIndex] == self.stdout.subdict.boundary:

                wstr  = "( "
                for vertice in self.faces[faceIndex]:
                    wstr += str(vertice) + " "
                wstr += ")"
                self.stdout.write(wstr)

        # ------------------------------------------------------------------- #

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

    # ----------------------------------------------------------------------- #

    def __init__(self):

        self.stdout = self.stdoutDict()
        self.vertices = vertices(self.stdout)
        self.blocks = blocks(self.stdout, self.vertices)
        self.boundaryFaces = boundaryFaces(self.stdout, self.vertices, self.blocks)

    # ----------------------------------------------------------------------- #

    def boundarySubDict(self, name, typename):

        self.stdout.boundarySubDict(name, typename)
        return True

    # ----------------------------------------------------------------------- #

    def subDict(self, name):

        self.stdout.subDict(name)
        return True

    # ----------------------------------------------------------------------- #

    def header(self, scale):

        self.stdout.header(scale)

    # ----------------------------------------------------------------------- #

    def footer(self):

        self.stdout.footer()

# --------------------------------------------------------------------------- #
# --- End of module --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

