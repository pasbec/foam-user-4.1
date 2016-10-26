#!/usr/bin/python
# October 2016
# Vladimir Galindo (v.galindo@hzdr.de)
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import numpy as np

from foamTools.ioInfo import fileGetPath, objectIndent, objectHeader, objectFooter

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

H0 = 0.029999
R0 = 0.029999

nr = 61
nz = 61

print "H0, R0 :", str(H0) + ",", str(R0)
print "nr, nz : ", nr, ",", nz

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Define short indented line with line break
def i(iL,cS,eS='\n'): return objectIndent(cS + eS,iLevel=iL)

# --------------------------------------------------------------------------- #
# --- Main program sequence ------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Mesh and discretisation
rl = np.linspace(0.0,R0,nr)
zl = np.linspace(0.0,H0,nz)

with open(fileGetPath('lorentzForce.samleDict'),'w') as cf:

    cf.write(objectHeader('sampleDict', 'dictionary'))

    cf.write(i(0, 'interpolationScheme cellPointFace;'))
    cf.write(i(0, 'setFormat           raw;'))
    cf.write(i(0, 'surfaceFormat       null;'))
    cf.write(i(0, 'sets '))
    cf.write(i(0, '('))
    cf.write(i(1, 'planexz'))
    cf.write(i(1, '{'))
    cf.write(i(2, 'type cloud;'))
    cf.write(i(2, 'axis    xyz;'))
    cf.write(i(2, 'points'))
    cf.write(i(2, '('))

    z0 = 0.03

    for ri in range(len(rl)):
        for zi in range(len(zl)):

            cf.write(i(3,'(' + str(rl[ri]) + ' 0.0 ' + str(z0 + zl[zi]) + ')'))

    cf.write(i(2, ');'))
    cf.write(i(1, '}'))
    cf.write(i(0, ');'))

    cf.write('\n')

    cf.write(i(0, 'surfaces '))
    cf.write(i(0, '('))
    cf.write(i(0, ');'))

    cf.write('\n')

    cf.write(i(0, 'fields '))
    cf.write(i(0, '('))
    cf.write(i(1, 'FL'))
    cf.write(i(0, ');'))

    cf.write('\n')

    cf.write(objectFooter())

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
