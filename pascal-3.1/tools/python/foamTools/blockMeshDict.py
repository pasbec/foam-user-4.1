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

    fileName = None
    indentLevel = 0

    # ----------------------------------------------------------------------- #

    def _indent(self, level, string, end):

        if not (type(level) == int \
            and type(string) == str \
            and type(end) == str):

            raise KeyError()

        return objectIndent(string + end, iLevel=level)

    # ----------------------------------------------------------------------- #

    def _write(self, string):

        if not type(string) == str:

            raise KeyError()

        if self.fileName == None:

            sys.stdout.write(string)

        else:

            with open(self.fileName,'a') as f:

                f.write(string)

    # ----------------------------------------------------------------------- #

    def line(self):
        self._write("\n")

    # ----------------------------------------------------------------------- #

    def write(self, string, ind=True, end="\n"):

        if not (type(string) == str \
            and type(ind) == bool \
            and type(end) == str):

            raise KeyError()

        level = self.indentLevel
        if not ind: level = 0

        wstr = self._indent(level,str(string), end)
        self._write(wstr)



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

        baseFaces = [[0,1],[2,3],[4,5]]

        baseEdges = [[0,1,2,3],[4,5,6,7],[8,9,10,11]]

        baseTags = {"x": 0, "y": 1, "z": 2}

        baseGetTags = {0: "x", 1: "y", 2: "z"}

        vertices = [i for i in range(8)]

        verticeBaseBase = [[None,0,[0,1],1,2,[0,2],[0,1,2],[1,2]],
                           [0,None,1,[0,1],[0,2],2,[1,2],[0,1,2]],
                           [[0,1],1,None,0,[0,1,2],[1,2],2,[0,2]],
                           [1,[0,1],0,None,[1,2],[0,1,2],[0,2],2],
                           [2,[0,2],[0,1,2],[1,2],None,0,[0,1],1],
                           [[0,2],2,[1,2],[0,1,2],0,None,1,[0,1]],
                           [[0,1,2],[1,2],2,[0,2],[0,1],1,None,0],
                           [[1,2],[0,1,2],[0,2],2,1,[0,1],0,None]]

        verticeBaseSign = [[None,True,[True,True],True,True,[True,True],[True,True,True],[True,True]],
                           [False,None,True,[False,True],[False,True],True,[True,True],[False,True,True]],
                           [[False,False],False,None,False,[False,False,True],[False,True],True,[False,True]],
                           [False,[True,False],True,None,[False,True],[True,False,True],[True,True],True],
                           [False,[True,False],[True,True,False],[True,False],None,True,[True,True],True],
                           [[False,False],False,[True,False],[False,True,False],False,None,True,[False,True]],
                           [[False,False,False],[False,False],False,[False,False],[False,False],False,None,False],
                           [[False,False],[True,False,False],[True,False],False,False,[True,False],True,None]]

        faces = [i for i in range(6)]

        faceBase = [0,0,1,1,2,2]

        faceOpposites = [1,0,3,2,5,4]

        faceVertices = [[0,3,7,4],[1,2,6,5],
                        [0,1,5,4],[2,3,7,6],
                        [0,1,2,3],[4,5,6,7]]

        faceEdges = [[4,11,7,8],[5,10,6,9],
                     [0,9,3,8],[1,10,2,11],
                     [0,5,1,4],[3,6,2,7]]

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

        edges = [i for i in range(12)]

        edgeBase = [0,0,0,0,1,1,1,1,2,2,2,2]

        edgeVertices = [[0,1],[3,2],[7,6],[4,5],
                        [0,3],[1,2],[5,6],[4,7],
                        [0,4],[1,5],[2,6],[3,7]]

        edgeFaces = [[2,4],[3,4],[3,5],[2,5],
                     [0,4],[1,4],[1,5],[0,5],
                     [0,2],[1,2],[1,3],[0,3]]

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

    def _getPathVerticeBase(self,vertices):

        bases = []
        signs = []

        for i in range(len(vertices)):

            if i < (len(vertices)-1):

                bases.append(self.topo.verticeBaseBase[vertices[i]][vertices[i+1]])
                signs.append(self.topo.verticeBaseSign[vertices[i]][vertices[i+1]])

            else:

                bases.append(self.topo.verticeBaseBase[vertices[len(vertices)-1]][vertices[0]])
                signs.append(self.topo.verticeBaseSign[vertices[len(vertices)-1]][vertices[0]])

        return bases, signs


    # ----------------------------------------------------------------------- #

    def _getSharedVertices(self, vertices, otherVertices):

        if not (type(vertices) == list \
            and type(otherVertices) == list):

            raise KeyError()

        return list(set(vertices) & set(otherVertices))

    # ----------------------------------------------------------------------- #

    def _getNeighbours(self, block, otherBlock):

        if not (type(block) == int \
            and type(otherBlock) == int):

            raise KeyError()

        # Get both sets of block vertices
        blockVertices = self.blockVertices[block]
        otherBlockVertices = self.blockVertices[otherBlock]

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

    def _setNeighbourData(self, block):

        for otherBlock, otherBlockLabel in enumerate(self.labels):

            if not block == otherBlock:

                neighbours, neighbourType, neighbourVertices = \
                    self._getNeighbours(block, otherBlock)

                if neighbourType == "face":

                    try:
                        self.neighbours[block].index(otherBlock)
                    except:
                        if self.neighbours[block] == None:
                            self.neighbours[block] = []
                        self.neighbours[block].append(otherBlock)
                    try:
                        self.neighbours[otherBlock].index(block)
                    except:
                        if self.neighbours[otherBlock] == None:
                            self.neighbours[otherBlock] = []
                        self.neighbours[otherBlock].append(block)

                    # Identify all face neighbours
                    blockVertices = self.blockVertices[block]
                    otherBlockVertices = self.blockVertices[otherBlock]

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
                            self.faceNeighbours[block][localFaceIndex] = otherBlock
                        if len(otherSharedFaceVertices) == 4:
                            self.faceNeighbours[otherBlock][localFaceIndex] = block

                    # Identify all edge neighbours
                    blockVertices = self.blockVertices[block]
                    otherBlockVertices = self.blockVertices[otherBlock]

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
                                self.edgeNeighbours[block][localEdgeIndex].index(otherBlock)
                            except:
                                if self.edgeNeighbours[block][localEdgeIndex] == None:
                                    self.edgeNeighbours[block][localEdgeIndex] = []
                                self.edgeNeighbours[block][localEdgeIndex].append(otherBlock)
                        if len(otherSharedEdgeVertices) == 2:
                            try:
                                self.edgeNeighbours[otherBlock][localEdgeIndex].index(block)
                            except:
                                if self.edgeNeighbours[otherBlock][localEdgeIndex] == None:
                                    self.edgeNeighbours[otherBlock][localEdgeIndex] = []
                                self.edgeNeighbours[otherBlock][localEdgeIndex].append(block)

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

            block = self.labels.index(label)

            self.blockVertices[block] = blockVertices
            self.blockVerticeLabels[block] = verticeLabels
            self.divider[block] = divider
            self.gradings[block] = grading
            self.zones[block] = zone
            self.neighbours[block] = None
            self.faceNeighbours[block] = [None for i in range(6)]
            self.edgeNeighbours[block] = [None for i in range(12)]

        # New block
        except:

            block = len(self.labels)
            self.labels.append(label)
            self.labelIndex[label] = block

            self.blockVertices.append(blockVertices)
            self.blockVerticeLabels.append(verticeLabels)
            self.divider.append(divider)
            self.gradings.append(grading)
            self.zones.append(zone)
            self.neighbours.append(None)
            self.faceNeighbours.append([None for i in range(6)])
            self.edgeNeighbours.append([None for i in range(12)])

        self._setNeighbourData(block)

    # ----------------------------------------------------------------------- #

# TODO [High]: Split "setDivider" into smaller pieces

    def setDivider(self, blockLabels, par1, par2=None):

        processedBlocks = [False for i in self.labels]

        def setProcessedBlock(block):

            processedBlocks[block] = True

        def getProcessedBlock(block):

            return processedBlocks[block]

        def syncBlocks(block, blockDivider, base):

            print ">> CURRENT", block
            print

            faces = self.topo.faces
            baseFaces = self.topo.baseFaces[base]
            crossBaseFaces = \
                [f for i, f in enumerate(faces) if i not in baseFaces]

            print "   block =", block
            print "   blockDivide =", blockDivide
            print "   base =", base
            print "   baseFaces =", baseFaces
            print "   crossBaseFaces =", crossBaseFaces
            print
            print "   >> WRITE"
            print

            # Write
            self.divider[block][base] = blockDivider[base]

            setProcessedBlock(block)

            # Sync
            for face in crossBaseFaces:

                # Get face vertices
                faceVertices = self.topo.faceVertices[face]
                faceVerticesBaseBase, faceVerticesBaseSign = \
                    self._getPathVerticeBase(faceVertices)
                faceVerticeLabels = \
                    [ self.blockVertices[block][v] for v in faceVertices]

                # Get next neighbour and exit if there is none
                nextBlock = self.faceNeighbours[block][face]

                print "   face =", face
                print "   faceVertices =", faceVertices
                print "   faceVerticesBaseBase =", faceVerticesBaseBase
                print "   faceVerticesBaseSign =", faceVerticesBaseSign
                print "   faceVerticeLabels =", faceVerticeLabels
                print "   nextBlock =", nextBlock
                if not nextBlock == None: print "   getProcessedBlock(nextBlock) =", getProcessedBlock(nextBlock)
                print

                if nextBlock == None or \
                    getProcessedBlock(nextBlock) == True:

                    print "   >> CONTINUE"
                    print
                    continue

                # Get next face index
                nextFace = self.faceNeighbours[nextBlock].index(block)

                # Get next face vertices
                nextFaceVertices = \
                    self.topo.faceVertices[nextFace]
                nextFaceVerticesTransformed = \
                    [self.blockVertices[nextBlock].index(l) for l in faceVerticeLabels]
                nextFaceVerticeLabels = \
                    [ self.blockVertices[nextBlock][v] for v in nextFaceVertices]
                nextFaceVerticesBaseBase, nextFaceVerticesBaseSign = \
                    self._getPathVerticeBase(nextFaceVerticesTransformed)

                # Calculate base transformation map
                nextBaseOrientation = [None for i in self.topo.base]

                nextBaseOrientation[self.topo.faceBase[face]] = \
                    self.topo.faceBase[nextFace]
                nextBaseOrientation[faceVerticesBaseBase[0]] = \
                    nextFaceVerticesBaseBase[0]
                nextBaseOrientation[faceVerticesBaseBase[1]] = \
                    nextFaceVerticesBaseBase[1]

                # Calculate face transformation map
                nextFaceOrientation = [None for i in self.topo.faces]

                nextFaceOrientation[self.topo.faceOpposites[face]] = nextFace
                nextFaceOrientation[face] = self.topo.faceOpposites[nextFace]
# TODO [High]: Currently we do no calculate nextFaceOrientation completely! Fix this

                # Assemble block divider for next block
                nextBlockDivider = [blockDivider[b] for b in nextBaseOrientation]

                # Extract neighbour divide
                nextBase = nextBaseOrientation[base]

                print "   nextFace =", nextFace
                print "   nextFaceVertices =", nextFaceVertices
                print "   nextFaceVerticesTransformed =", nextFaceVerticesTransformed
                print "   nextFaceVerticeLabels =", nextFaceVerticeLabels
                print "   nextFaceVerticesBaseBase =", nextFaceVerticesBaseBase
                print "   nextFaceVerticesBaseSign =", nextFaceVerticesBaseSign
                print "   nextBaseOrientation =", nextBaseOrientation
                print "   nextFaceOrientation =", nextFaceOrientation
                print "   nextBlockDivider =", nextBlockDivider
                print "   nextBase =", nextBase
                print

                print "   >> HOP"
                print

                syncBlocks(nextBlock, nextBlockDivider, nextBase)

        # ------------------------------------------------------------------- #

        if type(blockLabels) == int: blockLabels = [ blockLabels ]

        if not type(blockLabels) == list: raise KeyError()

        blocks = [ self.labelIndex[l] for l in blockLabels ]

        for block in blocks:

            # Divider given
            if not par2:

                if type(par1) == list:

                    blockDivider = par1

                    # Set divider for current block
                    self.divider[block] = blockDivider

                    # Sync blocks in all directions
                    for b in self.topo.base:
                        syncBlocks(block, self.divider[block], b)

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
                self.divider[block][blockDivideBase] = blockDivide

                # Sync blocks in direction of blockDivideBase
                syncBlocks(block, self.divider[block], blockDivideBase)

    # ----------------------------------------------------------------------- #

    def write(self):

        for block, blockLabel in enumerate(self.labels):

            self.stdout.write("hex", end=" ")

        # ------------------------------------------------------------------- #

            wstr  = "( "
            for vertice in self.blockVertices[block]:
                wstr += str(vertice) + " "
            wstr += ")"

            self.stdout.write(wstr, ind=False, end=" ")

        # ------------------------------------------------------------------- #

            if self.zones[block]:

                self.stdout.write(str(self.zones[block]), ind=False, end=" ")

        # ------------------------------------------------------------------- #

            wstr  = "( "
            for divide in self.divider[block]:
                wstr += str(divide) + " "
            wstr += ")"

            self.stdout.write(wstr, ind=False, end=" ")

        # ------------------------------------------------------------------- #

            if len(self.gradings[block]) == 3:

                self.stdout.write("simpleGrading", ind=False, end=" ")

            elif len(self.gradings[block]) == 12:

                self.stdout.write("edgeGrading", ind=False, end=" ")

        # ------------------------------------------------------------------- #

            wstr  = "( "
            for grading in self.gradings[block]:
                wstr += str(grading) + " "
            wstr += ")"

            self.stdout.write(wstr, ind=False)

    # ----------------------------------------------------------------------- #

    def printData(self, blockLabels=None):

        # Only print data for given block
        if blockLabels:

            if type(blockLabels) == int: blockLabels = [ blockLabels ]

            if not type(blockLabels) == list: raise KeyError()

            printBlocks = [ self.labelIndex[l] for l in blockLabels ]

        # Print all blocks
        else:

            printBlocks = [ self.labelIndex[l] for l in self.labels ]

        # ------------------------------------------------------------------- #

        for block in printBlocks:

            print "index:", block
            print "label:", self.labels[block]
            print "vertices: ", self.blockVertices[block]
            print "verticeLabels: ", self.blockVerticeLabels[block]
            print "divider: ", self.divider[block]
            print "gradings: ", self.gradings[block]
            print "zone: ", self.zones[block]
            print "neighbours: ", self.neighbours[block]
            print "faceNeighbours: ", self.faceNeighbours[block]
            print "edgeNeighbours: ", self.edgeNeighbours[block]
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

        def __init__(self, fileName=None):

            self.fileName = fileName

            if not fileName == None:

                with open(self.fileName,'w') as f: f.write("")

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
            self.write(wstr)
            wstr = "convertToMeters" + " " + str(scale) + ";"
            self.write(wstr)
            self.line()

        def footer(self):

            if self.subdict.opened: self.subDictFooter()
            wstr = objectFooter()
            self.write(wstr)

    # ----------------------------------------------------------------------- #

    def __init__(self, fileName=None):

        self.stdout = self.stdoutDict(fileName)
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

    # ----------------------------------------------------------------------- #

    def manual(self, string, end="\n"):

        self.stdout.write(string, end=end)

# --------------------------------------------------------------------------- #
# --- End of module --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

