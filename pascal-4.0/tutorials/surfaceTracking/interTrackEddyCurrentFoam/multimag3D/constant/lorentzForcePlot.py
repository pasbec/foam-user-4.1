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

fontsize   = 14
fontfamily = 'serif'

print "nr, nz     : ", nr, ",", nz
print "fontsize   : ", fontsize
print "fontfamily : ", fontfamily

setBaseName = 'lorentzForce'

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

data = dict()

R = dict()
Z = dict()

Fr = dict()
Fa = dict()
Fz = dict()

F  = dict()

# --------------------------------------------------------------------------- #
# --- Datasets -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

set = 'Analytical'
data[set] = np.genfromtxt(fileGetPath(setBaseName+set+'.dat'), comments='#')

R[set]  = data[set][:,0].reshape(nr,nz)
Z[set]  = data[set][:,1].reshape(nr,nz)
Fa[set] = data[set][:,2].reshape(nr,nz)

# Scale to mm
for ri in range(nr):
    for zi in range(nz):

        R[set][ri,zi] = 1000.0 * (R[set][ri,zi])

        Z[set][ri,zi] = 1000.0 * (Z[set][ri,zi])



set = 'Opera3D'
data[set] = np.genfromtxt(fileGetPath(setBaseName+set+'.dat'), comments='#')

R[set]    = data[set][:,0].reshape(nr,nz)
Z[set]    = data[set][:,1].reshape(nr,nz)
Fr[set]   = data[set][:,2].reshape(nr,nz)
Fa[set]   = data[set][:,3].reshape(nr,nz)
Fz[set]   = data[set][:,4].reshape(nr,nz)
F[set]    = np.zeros(R[set].shape)

# Scale to mm and extract force
for ri in range(nr):
    for zi in range(nz):

        F[set][ri,zi] = m.sqrt(
            Fr[set][ri,zi]**2 + Fa[set][ri,zi]**2 + Fz[set][ri,zi]**2)



set = 'EddyCurrentFoam'
data[set] = np.genfromtxt(fileGetPath(setBaseName+set+'.dat'), comments='#')

R[set]  = data[set][:,0].reshape(nr,nz)
Z[set]  = data[set][:,2].reshape(nr,nz)
Fr[set] = data[set][:,3].reshape(nr,nz)
Fa[set] = data[set][:,4].reshape(nr,nz)
Fz[set] = data[set][:,5].reshape(nr,nz)
F[set]  = np.zeros(R[set].shape)

z0 = 0.03

# Scale to mm, turn y-axis and extract force
for ri in range(nr):
    for zi in range(nz):

        R[set][ri,zi] = 1000.0 * (R[set][ri,zi])

        Z[set][ri,zi] = 1000.0 * (Z[set][ri,zi] - z0)

        Fa[set][ri,zi] = -Fa[set][ri,zi]

        F[set][ri,zi] = m.sqrt(
            Fr[set][ri,zi]**2 + Fa[set][ri,zi]**2 + Fz[set][ri,zi]**2)

# --------------------------------------------------------------------------- #
# --- Common plot settings -------------------------------------------------- #
# --------------------------------------------------------------------------- #

rLabel = '$r ~ [\mathrm{mm}]$'
zLabel = '$z ~ [\mathrm{mm}]$'

FrLevels = np.linspace(-0.30,-0.05,6)
FaLevels = np.linspace( 0.20, 2.00,10)
FzLevels = np.linspace(-0.40,-0.05,8)
Flevels  = FaLevels
colors = 'black'

latexify(fontsize=fontsize, fontfamily=fontfamily)

#fig = dict()

# --------------------------------------------------------------------------- #
# --- Plots ----------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

fig = plt.figure(0)
axe = fig.add_subplot(111)

axe.set_xlim([0,30])
axe.set_ylim([0,30])

axe.set_xlabel(rLabel)
axe.set_ylabel(zLabel)

axe.set_aspect('equal')

ct0 = axe.contour(
    R['Analytical'], Z['Analytical'], Fa['Analytical'],
    levels=FaLevels, colors=colors
    )

cl0 = axe.clabel(ct0, ct0.levels[0::2], inline=True, fmt='%g', fontsize=fontsize)

fig.savefig(fileGetPath('lorentzForceAnalytical.pdf'),bbox_inches="tight")



fig = plt.figure(1)
axe = fig.add_subplot(111)

axe.set_xlim([0,30])
axe.set_ylim([0,30])

axe.set_xlabel(rLabel)
axe.set_ylabel(zLabel)

axe.set_aspect('equal')

ct1 = axe.contour(
    R['Opera3D'], Z['Opera3D'], Fr['Opera3D'],
    levels=FrLevels, colors=colors, linestyles='dashed'
    )

ct2 = axe.contour(
    R['EddyCurrentFoam'], Z['EddyCurrentFoam'], Fr['EddyCurrentFoam'],
    levels=FrLevels, colors=colors, linestyles='solid'
    )

cl2 = axe.clabel(ct2, ct2.levels[0::2], inline=True, fmt='%g', fontsize=fontsize)
[txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=2)) for txt in cl2]

fig.savefig(fileGetPath('lorentzForceComparisonFr.pdf'),bbox_inches="tight")



fig = plt.figure(2)
axe = fig.add_subplot(111)

axe.set_xlim([0,30])
axe.set_ylim([0,30])

axe.set_xlabel(rLabel)
axe.set_ylabel(zLabel)

axe.set_aspect('equal')

ct0 = axe.contour(
    R['Analytical'], Z['Analytical'], Fa['Analytical'],
    levels=FaLevels, colors=colors, linestyles='dotted'
    )

ct1 = axe.contour(
    R['Opera3D'], Z['Opera3D'], Fa['Opera3D'],
    levels=FaLevels, colors=colors, linestyles='dashed'
    )

ct2 = axe.contour(
    R['EddyCurrentFoam'], Z['EddyCurrentFoam'], Fa['EddyCurrentFoam'],
    levels=FaLevels, colors=colors, linestyles='solid'
    )

cl2 = axe.clabel(ct2, ct2.levels[0::2], inline=True, fmt='%g', fontsize=fontsize)
[txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=2)) for txt in cl2]

fig.savefig(fileGetPath('lorentzForceComparisonFa.pdf'),bbox_inches="tight")



fig = plt.figure(3)
axe = fig.add_subplot(111)

axe.set_xlim([0,30])
axe.set_ylim([0,30])

axe.set_xlabel(rLabel)
axe.set_ylabel(zLabel)

axe.set_aspect('equal')

ct1 = axe.contour(
    R['Opera3D'], Z['Opera3D'], Fz['Opera3D'],
    levels=FzLevels, colors=colors, linestyles='dashed'
    )

ct2 = axe.contour(
    R['EddyCurrentFoam'], Z['EddyCurrentFoam'], Fz['EddyCurrentFoam'],
    levels=FzLevels, colors=colors, linestyles='solid'
    )

cl2 = axe.clabel(ct2, ct2.levels[0::2], inline=True, fmt='%g', fontsize=fontsize)
[txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=2)) for txt in cl2]

fig.savefig(fileGetPath('lorentzForceComparisonFz.pdf'),bbox_inches="tight")



fig = plt.figure(4)
axe = fig.add_subplot(111)

axe.set_xlim([0,30])
axe.set_ylim([0,30])

axe.set_xlabel(rLabel)
axe.set_ylabel(zLabel)

axe.set_aspect('equal')

ct1 = axe.contour(
    R['Opera3D'], Z['Opera3D'], F['Opera3D'],
    levels=Flevels, colors=colors, linestyles='dashed'
    )

ct2 = axe.contour(
    R['EddyCurrentFoam'], Z['EddyCurrentFoam'], F['EddyCurrentFoam'],
    levels=Flevels, colors=colors, linestyles='solid'
    )

cl2 = axe.clabel(ct2, ct2.levels[0::2], inline=True, fmt='%g', fontsize=fontsize)
[txt.set_bbox(dict(facecolor='white', edgecolor='none', pad=2)) for txt in cl2]

fig.savefig(fileGetPath('lorentzForceComparisonF.pdf'),bbox_inches="tight")

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
