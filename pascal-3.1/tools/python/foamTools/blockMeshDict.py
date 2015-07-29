#!/usr/bin/python
# Module with helper functions to write a blockdict with python
# July 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

from foamTools.ioInfo import *

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

subDictOpen = False
boundarySubDictOpen = False
indentLevel = 0



def printi(s, ind=True, end="\n"):
    def i(iL,cS,eS='\n'): return objectIndent(cS + eS,iLevel=iL)
    iL = indentLevel
    if not ind: iL = 0
    sys.stdout.write(i(iL,str(s), end))



def pLine():
    print ""



def pBoundarySubDictHeader(n,t):
    global boundarySubDictOpen, indentLevel
    boundarySubDictOpen = True
    printi(str(n))
    printi("{")
    indentLevel += 1
    printi("type" + " " + str(t) + ";")
    printi("faces")
    printi("(")
    indentLevel += 1



def pBoundarySubDictFooter(l = True):
    global boundarySubDictOpen, indentLevel
    boundarySubDictOpen = False
    indentLevel -= 1
    printi(");")
    indentLevel -= 1
    printi("}")
    if l: pLine()



def pBoundarySubDict(n,t):
    global subDictOpen, boundarySubDictOpen
    if boundarySubDictOpen: pBoundarySubDictFooter()
    pBoundarySubDictHeader(n,t)
    return True



def pSubDictHeader(t):
    global subDictOpen, indentLevel
    subDictOpen = True
    printi(str(t))
    printi("(")
    indentLevel += 1



def pSubDictFooter(l = True):
    global subDictOpen, indentLevel
    subDictOpen = False
    indentLevel -= 1
    printi(");")
    if l: pLine()



def pSubDict(t):
    global boundarySubDictOpen, subDictOpen
    if boundarySubDictOpen: pBoundarySubDictFooter(False)
    if subDictOpen: pSubDictFooter()
    pSubDictHeader(t)
    return True



def pHeader(s):
    sys.stdout.write(objectHeader("blockMeshDict", "dictionary"))
    ss = "convertToMeters" + " " + str(s) + ";"
    printi(ss)
    pLine()



def pFooter():
    global subDictOpen
    if subDictOpen: pSubDictFooter()
    sys.stdout.write(objectFooter())



def pVert(i,p):

    ps  = "(" + str(float(p[0])) + " "
    ps +=       str(float(p[1])) + " "
    ps +=       str(float(p[2])) + ")"
    printi(ps)



def pBlock(i,v,n,g=[1,1,1],z=""):

    def prints(s): printi(s,ind=False,end="")

    printi("hex", end="")

    prints(" ")

    vs  = "(" + str(v[0]) + " " + str(v[1]) + " "
    vs +=       str(v[2]) + " " + str(v[3]) + " "
    vs +=       str(v[4]) + " " + str(v[5]) + " "
    vs +=       str(v[6]) + " " + str(v[7]) + ")"
    prints(vs)

    prints(" ")

    if z: prints(z + " ")

    ns = "(" + str(n[0]) + " " + str(n[1]) + " " + str(n[2]) + ")"
    prints(ns)

    prints(" ")

    prints("simpleGrading")

    prints(" ")

    gs = "(" + str(g[0]) + " " + str(g[1]) + " " + str(g[2]) + ")"
    printi(gs, ind=False)




def pFace(f):

    fs  = "(" + str(f[0]) + " " + str(f[1]) + " "
    fs +=       str(f[2]) + " " + str(f[3]) + ")"
    printi(fs)

# --------------------------------------------------------------------------- #
# --- End of module --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

