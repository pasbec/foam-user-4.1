#!/usr/bin/python
# October 2016
# Vladimir Galindo (v.galindo@hzdr.de)
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import math as m
import numpy as np

import matplotlib.pyplot as plt

from foamTools.ioInfo import fileGetPath, objectIndent, objectHeader, objectFooter

from plotTools.latex import latexify

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

nr = 61
nz = 61

fontsize=12
fontfamily='serif'

print "nr, nz     : ", nr, ",", nz
print "fontsize   : ", fontsize
print "fontfamily : ", fontfamily

# --------------------------------------------------------------------------- #
# --- Main program sequence ------------------------------------------------- #
# --------------------------------------------------------------------------- #

latexify(fontsize=fontsize, fontfamily=fontfamily)



data0 = np.genfromtxt(fileGetPath('lorentzForceAnalytical.dat'), comments='#')

R0 = data0[:,0].reshape(nr,nz)
Z0 = data0[:,1].reshape(nr,nz)
F0 = data0[:,2].reshape(nr,nz)

# Scale to mm
for ri in range(nr):
    for zi in range(nz):

        R0[ri,zi] = 1000.0 * (R0[ri,zi])

        Z0[ri,zi] = 1000.0 * (Z0[ri,zi])

fig0, ax0 = plt.subplots()

ax0.set_xlim([0,30])
ax0.set_ylim([0,30])

ax0.set_aspect('equal')

#fc0 = ax0.contour(R0, Z0, F0, colors='black', levels=np.linspace(0.2,2.0,10))
fc0 = ax0.contour(R0, Z0, F0, levels=np.linspace(0.2,2.0,10))

ax0.clabel(fc0, fc0.levels[0::2], inline=True, fmt='%g', fontsize=fontsize)

fig0.savefig(fileGetPath('lorentzForceAnalytical.pdf'),bbox_inches="tight")



data1 = np.genfromtxt(fileGetPath('lorentzForce.dat'), comments='#')

R1  = data1[:,0].reshape(nr,nz)
Z1  = data1[:,2].reshape(nr,nz)
F1x = data1[:,3].reshape(nr,nz)
F1y = data1[:,4].reshape(nr,nz)
F1z = data1[:,5].reshape(nr,nz)
F1  = np.zeros(R1.shape)

z0 = 0.03

# Scale to mm and extract force
for ri in range(nr):
    for zi in range(nz):

        R1[ri,zi] = 1000.0 * (R1[ri,zi])

        Z1[ri,zi] = 1000.0 * (Z1[ri,zi] - z0)

        # Magnitude
        #F1[ri,zi] = m.sqrt(F1x[ri,zi]**2 + F1y[ri,zi]**2 + F1z[ri,zi]**2)

        # Phi-component
        F1[ri,zi] = abs(F1y[ri,zi])

fig1, ax1 = plt.subplots()

ax1.set_xlim([0,30])
ax1.set_ylim([0,30])

ax1.set_aspect('equal')

#fc1 = ax1.contour(R1, Z1, F1, colors='black', levels=np.linspace(0.2,2.0,10))
fc1 = ax1.contour(R1, Z1, F1, levels=np.linspace(0.2,2.0,10))

ax1.clabel(fc1, fc1.levels[0::2], inline=True, fmt='%g', fontsize=fontsize)

fig1.savefig(fileGetPath('lorentzForce.pdf'),bbox_inches="tight")


#plt.show(block=False)
#time.sleep(10)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
