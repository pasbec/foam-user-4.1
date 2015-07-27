#!/usr/bin/python
import math as m

## * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * ##

subDictOpen = None
patchSubDictOpen = None

def pLine():
    print ""

def pBoundarySubDictHeader(n,t):
    print "    %s" % (n)
    print "    {"
    print "        type %s;" % (t)
    print "        faces"
    print "        ("

def pBoundarySubDictFooter(l = True):
    print "        );"
    print "    }"
    if l:
        pLine()

def pBoundarySubDict(n,t):
    global patchSubDictOpen
    if patchSubDictOpen:
        pBoundarySubDictFooter()
    pBoundarySubDictHeader(n,t)
    patchSubDictOpen = True

def pSubDictHeader(t):
    print "%s" % (t)
    print "("

def pSubDictFooter(l = True):
    print ");"
    if l:
        pLine()

def pSubDict(t):
    global patchSubDictOpen
    if patchSubDictOpen:
        pBoundarySubDictFooter(False)
    global subDictOpen
    if subDictOpen:
        pSubDictFooter()
    pSubDictHeader(t)
    subDictOpen = True

def pHeader():
    global subDictOpen
    subDictOpen = False
    print """\
/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | foam-extend: Open Source CFD                    |
|  \\\\    /   O peration     | Version:     3.1                                |
|   \\\\  /    A nd           | Web:         http://www.extend-project.de       |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version 2.0;
    format  ascii;
    class   dictionary;
    object  blockMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\
"""
    pLine()

def pFooter():
    global subDictOpen
    if subDictOpen:
        pSubDictFooter()
    print """\
// ************************************************************************* //
"""

def pConvertToMeters(s):
    print "convertToMeters %g;" % (s)
    pLine()

def pVert(n,x,y,z):
    print "    (%g %g %g)" % (x, y, z)

## ************************************************************************* ##
