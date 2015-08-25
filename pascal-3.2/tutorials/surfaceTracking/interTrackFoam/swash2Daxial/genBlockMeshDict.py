#!/usr/bin/python
import math as m

nScale = 0.5

x1   = 85.0
x11  = x1/2.0
x12  = x1/2.0

nx1  = m.ceil(nScale*x1 )
nx11 = m.ceil(nScale*x11)
nx12 = m.ceil(nScale*x12)

z0   = 0.0
z1   = 235.0-1.0/(1.0*nScale)
z11  = z1/2.0
z12  = z1/2.0
z2   = 235.0-z1

nz1  = m.ceil(nScale*z1 )
nz11 = m.ceil(nScale*z11)
nz12 = m.ceil(nScale*z12)
nz2  = m.ceil(nScale*z2 )

phiGrad  = 5.0
phi      = m.pi * phiGrad/180.0
phiHalf  = phi/2.0
pSin     = m.sin(phiHalf)
pCos     = m.cos(phiHalf)

def pvert(x,y,z,n):
  print "  (%f %f %f) \t\t// vertex %d" % (x, y, z, n)

print """//
FoamFile
{
  version 2.0;
  format  ascii;
  class   dictionary;
  object  blockMeshDict;
}

convertToMeters 0.001;

vertices
("""

pvert(0       , 0       , z0   , 0)
pvert(pCos*x11,-pSin*x11, z0   , 1)
pvert(pCos*x1 ,-pSin*x1 , z0   , 2)

pvert(0       , 0       , z11  , 3)
pvert(pCos*x11,-pSin*x11, z11  , 4)
pvert(pCos*x1 ,-pSin*x1 , z11  , 5)

pvert(0       , 0       , z1   , 6)
pvert(pCos*x11,-pSin*x11, z1   , 7)
pvert(pCos*x1 ,-pSin*x1 , z1   , 8)

pvert(pCos*x11, pSin*x11, z0   , 9)
pvert(pCos*x1 , pSin*x1 , z0   ,10)

pvert(pCos*x11, pSin*x11, z11  ,11)
pvert(pCos*x1 , pSin*x1 , z11  ,12)

pvert(pCos*x11, pSin*x11, z1   ,13)
pvert(pCos*x1 , pSin*x1 , z1   ,14)

pvert(0       , 0       , z1+z2,15)
pvert(pCos*x11,-pSin*x11, z1+z2,16)
pvert(pCos*x1 ,-pSin*x1 , z1+z2,17)

pvert(pCos*x11, pSin*x11, z1+z2,18)
pvert(pCos*x1 , pSin*x1 , z1+z2,19)

print """);

blocks
("""

print "  hex (0 1 9 0 3 4 11 3)      liCond (%d %d %d) simpleGrading (1 1 1)" % (m.ceil(nx11), 1, m.ceil(nz11)) # block 0
print "  hex (1 2 10 9 4 5 12 11)    liCond (%d %d %d) simpleGrading (1 1 1)" % (m.ceil(nx12), 1, m.ceil(nz11)) # block 1
print "  hex (3 4 11 3 6 7 13 6)     liCond (%d %d %d) simpleGrading (1 1 1)" % (m.ceil(nx11), 1, m.ceil(nz12)) # block 2
print "  hex (4 5 12 11 7 8 14 13)   liCond (%d %d %d) simpleGrading (1 1 1)" % (m.ceil(nx12), 1, m.ceil(nz12)) # block 3
print "  hex (6 7 13 6 15 16 18 15)  liCond (%d %d 1 ) simpleGrading (1 1 1)" % (m.ceil(nx11), 1              ) # block 4
print "  hex (7 8 14 13 16 17 19 18) liCond (%d %d 1 ) simpleGrading (1 1 1)" % (m.ceil(nx12), 1              ) # block 5


print """);

edges
("""

print """);

boundary
(
  freeSurface
  {
    type patch;
    faces
    (
      (15 16 18 15)
      (16 17 19 18)
    );
  }
  side
  {
    type wall;
    faces
    (
      ( 8 14 19 17)
      ( 2 10 12  5)
      ( 5 12 14  8)
    );
  }
  bottom
  {
    type wall;
    faces
    (
      ( 0  1  9  0)
      ( 1  2 10  9)
    );
  }
  front
  {
    type wedge;
    faces
    (
      ( 0  1  4  3)
      ( 1  2  5  4)
      ( 3  4  7  6)
      ( 4  5  8  7)
      ( 6  7 16 15)
      ( 7  8 17 16)
    );
  }
  back
  {
    type wedge;
    faces
    (
      ( 0  3 11  9)
      ( 9 11 12 10)
      ( 3  6 13 11)
      (11 13 14 12)
      ( 6 13 18 15)
      (13 14 19 18)
    );
  }
  axis
  {
    type empty;
    faces
    (
      ( 0  3  3  0)
      ( 3  6  6  3)
    );
  }
);

mergePatchPairs
(
);"""
