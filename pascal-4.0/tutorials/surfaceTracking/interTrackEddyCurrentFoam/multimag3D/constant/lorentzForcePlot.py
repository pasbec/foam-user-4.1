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

from foamTools.ioInfo import fileGetPath

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

baseName = 'lorentzForce'

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

data = dict()

R  = dict()
Z  = dict()

Fr = dict()
Fa = dict()
Fz = dict()
F  = dict()

FaE = dict()
FzE = dict()
FrE = dict()
FE  = dict()

# --------------------------------------------------------------------------- #

set = 'Analytical'

data[set] = np.genfromtxt(fileGetPath(baseName+set+'.dat'), comments='#')

R[set]  = data[set][:,0].reshape(nr,nz)
Z[set]  = data[set][:,1].reshape(nr,nz)

Fa[set] = data[set][:,2].reshape(nr,nz)

# Scale to mm
for ri in range(nr):
    for zi in range(nz):

        R[set][ri,zi] = 1000.0 * (R[set][ri,zi])

        Z[set][ri,zi] = 1000.0 * (Z[set][ri,zi])



set = 'Opera3D'

data[set] = np.genfromtxt(fileGetPath(baseName+set+'.dat'), comments='#')

R[set]  = data[set][:,0].reshape(nr,nz)
Z[set]  = data[set][:,1].reshape(nr,nz)

Fr[set] = data[set][:,2].reshape(nr,nz)
Fa[set] = data[set][:,3].reshape(nr,nz)
Fz[set] = data[set][:,4].reshape(nr,nz)
F[set]  = np.zeros(R[set].shape)

# Scale to mm and extract force
for ri in range(nr):
    for zi in range(nz):

        F[set][ri,zi] = m.sqrt(
            Fr[set][ri,zi]**2 + Fa[set][ri,zi]**2 + Fz[set][ri,zi]**2)



set = 'EddyCurrentFoam'

data[set] = dict()

R[set]  = dict()
Z[set]  = dict()

Fr[set] = dict()
Fa[set] = dict()
Fz[set] = dict()

F[set]  = dict()

meshes = ['0.125', '0.250', '0.375', '0.500', '0.750', '1.000', '1.500', '2.000']
for mesh in meshes:

    data[set][mesh] = np.genfromtxt(fileGetPath(baseName+set+'_'+mesh+'.dat'), comments='#')

    d = data[set][mesh]

    R[set][mesh]  = data[set][mesh][:,0].reshape(nr,nz)
    Z[set][mesh]  = data[set][mesh][:,2].reshape(nr,nz)

    Fr[set][mesh] = data[set][mesh][:,3].reshape(nr,nz)
    Fa[set][mesh] = data[set][mesh][:,4].reshape(nr,nz)
    Fz[set][mesh] = data[set][mesh][:,5].reshape(nr,nz)
    F[set][mesh]  = np.zeros(R[set][mesh].shape)

    # Scale to mm, flip y-axis and extract force
    for ri in range(nr):
        for zi in range(nz):

            R[set][mesh][ri,zi] = 1000.0 * (R[set][mesh][ri,zi])

            Z[set][mesh][ri,zi] = 1000.0 * (Z[set][mesh][ri,zi] - 0.03)

            Fa[set][mesh][ri,zi] = -Fa[set][mesh][ri,zi]

            F[set][mesh][ri,zi] = m.sqrt(
                Fr[set][mesh][ri,zi]**2 + Fa[set][mesh][ri,zi]**2 + Fz[set][mesh][ri,zi]**2)

FrE[set] = dict()
FaE[set] = dict()
FzE[set] = dict()

FE[set]  = dict()

for smesh in meshes[:-1]:

    FrE[set][smesh]  = np.zeros(R[set][mesh].shape)
    FaE[set][smesh]  = np.zeros(R[set][mesh].shape)
    FzE[set][smesh]  = np.zeros(R[set][mesh].shape)
    FE[set][smesh]  = np.zeros(R[set][mesh].shape)

    # Scale to mm, flip y-axis and extract force
    for ri in range(nr):
        for zi in range(nz):

            FrE[set][smesh][ri,zi] = abs(Fr[set][meshes[-1]][ri,zi] - Fr[set][mesh][ri,zi])
            FaE[set][smesh][ri,zi] = abs(Fa[set][meshes[-1]][ri,zi] - Fa[set][mesh][ri,zi])
            FzE[set][smesh][ri,zi] = abs(Fz[set][meshes[-1]][ri,zi] - Fz[set][mesh][ri,zi])
            FE[set][smesh][ri,zi]  = abs(F[set][meshes[-1]][ri,zi] - F[set][mesh][ri,zi])

# --------------------------------------------------------------------------- #
# --- Error plot ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #



# --------------------------------------------------------------------------- #
# --- Contour plots --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

labelR = '$r ~ [\mathrm{mm}]$'
labelZ = '$z ~ [\mathrm{mm}]$'

levelsFr = np.linspace(-0.30,-0.05,6)
levelsFa = np.linspace( 0.20, 2.00,10)
levelsFz = np.linspace(-0.40,-0.05,8)
levelsF  = levelsFa

colors = 'black'

latexify(fontsize=fontsize, fontfamily=fontfamily)

plots = dict()

# --------------------------------------------------------------------------- #

def fig(p, name):

    p[name] = {'fig': plt.figure(), 'axes': dict()}
    f = p[name]

    fig = f['fig']
    axes = f['axes']

    def axe(f, axes, name):

        axes[name] = {'axe': fig.add_subplot(111), 'elements': dict()}
        a = axes[name]

        axe = a['axe']
        elem = a['elements']

        axe.set_xlim([0,30])
        axe.set_ylim([0,30])

        axe.set_xlabel(labelR)
        axe.set_ylabel(labelZ)

        axe.set_aspect('equal')

        def ele(f, a, elements, name):

            elements[name] = a.contour(
                R['Analytical'], Z['Analytical'], Fa['Analytical'],
                levels=levelsFa, colors=colors)

            e = elements[name]

            a.clabel(
                e, e.levels[0::2],
                inline=True, fmt='%g', fontsize=fontsize)

        ele(fig, axe, elem, 'contour')

    axe(fig, axes, 'Fa')

    fig.savefig(fileGetPath(baseName+name+'.pdf'), bbox_inches="tight")

fig(plots, 'Analytical')



def fig(p, name):

    p[name] = {'fig': plt.figure(), 'axes': dict()}
    f = p[name]

    fig = f['fig']
    axes = f['axes']

    def axe(f, axes, name):

        axes[name] = {'axe': fig.add_subplot(111), 'elements': dict()}
        a = axes[name]

        axe = a['axe']
        elem = a['elements']

        axe.set_xlim([0,30])
        axe.set_ylim([0,30])

        axe.set_xlabel(labelR)
        axe.set_ylabel(labelZ)

        axe.set_aspect('equal')

        def ele(f, a, elements, name):

            elements[name] = a.contour(
                R['Opera3D'], Z['Opera3D'], Fr['Opera3D'],
                levels=levelsFr, colors=colors, linestyles='dashed')

        ele(fig, axe, elem, 'cOpera3D')

        def ele(f, a, elements, name):

            set = 'EddyCurrentFoam'
            mesh = '1.000'

            elements[name] = a.contour(
                R[set][mesh], Z[set][mesh], Fr[set][mesh],
                levels=levelsFr, colors=colors, linestyles='solid')

            e = elements[name]

            cl = a.clabel(
                e, e.levels[0::2],
                inline=True, fmt='%g', fontsize=fontsize)

            [l.set_bbox(dict(facecolor='white', edgecolor='none', pad=2)) for l in cl]

        ele(fig, axe, elem, 'cEddyCurrentFoam')

    axe(fig, axes,'Fr')

    fig.savefig(fileGetPath(baseName+name+'.pdf'), bbox_inches="tight")

fig(plots, 'ComparisonFr')



def fig(p, name):

    p[name] = {'fig': plt.figure(), 'axes': dict()}
    f = p[name]

    fig = f['fig']
    axes = f['axes']

    def axe(f, axes, name):

        axes[name] = {'axe': fig.add_subplot(111), 'elements': dict()}
        a = axes[name]

        axe = a['axe']
        elem = a['elements']

        axe.set_xlim([0,30])
        axe.set_ylim([0,30])

        axe.set_xlabel(labelR)
        axe.set_ylabel(labelZ)

        axe.set_aspect('equal')

        def ele(f, a, elements, name):

            elements[name] = a.contour(
                R['Analytical'], Z['Analytical'], Fa['Analytical'],
                levels=levelsFa, colors=colors, linestyles='dotted')

        ele(fig, axe, elem, 'cAnalytical')

        def ele(f, a, elements, name):

            elements[name] = a.contour(
                R['Opera3D'], Z['Opera3D'], Fa['Opera3D'],
                levels=levelsFa, colors=colors, linestyles='dashed')

        ele(fig, axe, elem, 'cOpera3D')

        def ele(f, a, elements, name):

            set = 'EddyCurrentFoam'
            mesh = '1.000'

            elements[name] = a.contour(
                R[set][mesh], Z[set][mesh], Fa[set][mesh],
                levels=levelsFa, colors=colors, linestyles='solid')

            e = elements[name]

            cl = a.clabel(
                e, e.levels[0::2],
                inline=True, fmt='%g', fontsize=fontsize)

            [l.set_bbox(dict(facecolor='white', edgecolor='none', pad=2)) for l in cl]

        ele(fig, axe, elem, 'cEddyCurrentFoam')

    axe(fig, axes, 'Fa')

    fig.savefig(fileGetPath(baseName+name+'.pdf'), bbox_inches="tight")

fig(plots, 'ComparisonFa')



def fig(p, name):

    p[name] = {'fig': plt.figure(), 'axes': dict()}
    f = p[name]

    fig = f['fig']
    axes = f['axes']

    def axe(f, axes, name):

        axes[name] = {'axe': fig.add_subplot(111), 'elements': dict()}
        a = axes[name]

        axe = a['axe']
        elem = a['elements']

        axe.set_xlim([0,30])
        axe.set_ylim([0,30])

        axe.set_xlabel(labelR)
        axe.set_ylabel(labelZ)

        axe.set_aspect('equal')

        def ele(f, a, elements, name):

            elements[name] = a.contour(
                R['Opera3D'], Z['Opera3D'], Fz['Opera3D'],
                levels=levelsFz, colors=colors, linestyles='dashed')

        ele(fig, axe, elem, 'cOpera3D')

        def ele(f, a, elements, name):

            set = 'EddyCurrentFoam'
            mesh = '1.000'

            elements[name] = a.contour(
                R[set][mesh], Z[set][mesh], Fz[set][mesh],
                levels=levelsFz, colors=colors, linestyles='solid')

            e = elements[name]

            cl = a.clabel(
                e, e.levels[0::2],
                inline=True, fmt='%g', fontsize=fontsize)

            [l.set_bbox(dict(facecolor='white', edgecolor='none', pad=2)) for l in cl]

        ele(fig, axe, elem, 'cEddyCurrentFoam')

    axe(fig, axes, 'Fz')

    fig.savefig(fileGetPath(baseName+name+'.pdf'), bbox_inches="tight")

fig(plots, 'ComparisonFz')



def fig(p, name):

    p[name] = {'fig': plt.figure(), 'axes': dict()}
    f = p[name]

    fig = f['fig']
    axes = f['axes']

    def axe(f, axes, name):

        axes[name] = {'axe': fig.add_subplot(111), 'elements': dict()}
        a = axes[name]

        axe = a['axe']
        elem = a['elements']

        axe.set_xlim([0,30])
        axe.set_ylim([0,30])

        axe.set_xlabel(labelR)
        axe.set_ylabel(labelZ)

        axe.set_aspect('equal')

        def ele(f, a, elements, name):

            elements[name] = a.contour(
                R['Opera3D'], Z['Opera3D'], F['Opera3D'],
                levels=levelsF, colors=colors, linestyles='dashed')

        ele(fig, axe, elem, 'cOpera3D')

        def ele(f, a, elements, name):

            set = 'EddyCurrentFoam'
            mesh = '1.000'

            elements[name] = a.contour(
                R[set][mesh], Z[set][mesh], F[set][mesh],
                levels=levelsF, colors=colors, linestyles='solid')

            e = elements[name]

            cl = a.clabel(
                e, e.levels[0::2],
                inline=True, fmt='%g', fontsize=fontsize)

            [l.set_bbox(dict(facecolor='white', edgecolor='none', pad=2)) for l in cl]

        ele(fig, axe, elem, 'cEddyCurrentFoam')

    axe(fig, axes, 'F')

    fig.savefig(fileGetPath(baseName+name+'.pdf'), bbox_inches="tight")

fig(plots, 'ComparisonF')



#def fig(p, name):

    #p[name] = {'fig': plt.figure(), 'axes': dict()}
    #f = p[name]

    #fig = f['fig']
    #axes = f['axes']

    #axes['Fr'] = {'axe': fig.add_subplot(221), 'elements': p['ComparisonFr']['axes']['Fr']['elements']}
    #axes['Fa'] = {'axe': fig.add_subplot(222), 'elements': p['ComparisonFa']['axes']['Fa']['elements']}
    #axes['Fz'] = {'axe': fig.add_subplot(223), 'elements': p['ComparisonFz']['axes']['Fz']['elements']}
    #axes['F'] = {'axe': fig.add_subplot(224), 'elements': p['ComparisonF']['axes']['F']['elements']}

    #fig.savefig(fileGetPath(baseName+name+'.pdf'), bbox_inches="tight")

#fig(plots, 'Comparison')

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
