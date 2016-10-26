#!/usr/bin/python
# October 2016
# Vladimir Galindo (v.galindo@hzdr.de)
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import math
import numpy as np
import scipy as sp

import matplotlib as mpl
import matplotlib.pyplot as plt

from plotTools.latex import latexify

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

nr = 61
nz = 61

fontsize=14
fontfamily='serif'

print "nr, nz     : ", nr, ",", nz
print "fontsize   : ", fontsize
print "fontfamily : ", fontfamily

# --------------------------------------------------------------------------- #
# --- Main program sequence ------------------------------------------------- #
# --------------------------------------------------------------------------- #

data = np.genfromtxt('LorentzForce.dat', comments='#')

R = data[:,0].reshape(nr,nz)
Z = data[:,1].reshape(nr,nz)
F = data[:,2].reshape(nr,nz)

fig0, ax0 = plt.subplots()

ax0.set_xlim([0.0,0.03])
ax0.set_ylim([0.0,0.03])

ax0.set_aspect('equal')

#fc = ax0.contour(R, Z, F, colors='black', levels=np.linspace(0.2,2.0,10))
fc = ax0.contour(R, Z, F, levels=np.linspace(0.2,2.0,10))

ax0.clabel(fc, fc.levels[0::2], inline=True, fmt='%g', fontsize=fontsize)

latexify(fontsize=fontsize, fontfamily=fontfamily)

plt.savefig('LorentzForce.pdf',bbox_inches="tight")


#plt.show(block=False)
#time.sleep(10)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
