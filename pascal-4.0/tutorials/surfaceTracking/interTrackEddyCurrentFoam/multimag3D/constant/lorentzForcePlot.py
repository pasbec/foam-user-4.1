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
# --- Functions ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def norm_inf(r, e): return np.max(np.absolute(e)) / np.max(np.absolute(r))
def norm_1(r, e): return np.sum(np.absolute(e)) / np.sum(np.absolute(r))
def norm_2(r, e): return (np.sum(e**2) / np.sum(r**2))**0.5

norms = {'inf': norm_inf, '1': norm_1, '2': norm_2}

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

data = dict()

R = dict() # Radius
Z = dict() # Height

F = dict() # Force

E = dict() # Error
D = dict() # Delta of grid
N = dict() # Norm of error

# --------------------------------------------------------------------------- #

set = 'Analytical'

# Read data
if True:

    data[set] = np.genfromtxt(fileGetPath(baseName+set+'.dat'), comments='#')

    R[set]    = data[set][:,0].reshape(nr,nz)
    Z[set]    = data[set][:,1].reshape(nr,nz)

    F[set]    = [np.zeros(R[set].shape) for i in range(3)]
    F[set][1] = data[set][:,2].reshape(nr,nz)

    # Scale to mm
    R[set] = 1000.0 * R[set]
    Z[set] = 1000.0 * Z[set]

# --------------------------------------------------------------------------- #

set = 'Opera3D'

# Read data
if True:

    data[set] = np.genfromtxt(fileGetPath(baseName+set+'.dat'), comments='#')

    R[set]    = data[set][:,0].reshape(nr,nz)
    Z[set]    = data[set][:,1].reshape(nr,nz)

    F[set]    = [np.zeros(R[set].shape) for i in range(3)]
    F[set][0] = data[set][:,2].reshape(nr,nz)
    F[set][1] = data[set][:,3].reshape(nr,nz)
    F[set][2] = data[set][:,4].reshape(nr,nz)

# --------------------------------------------------------------------------- #

set = 'EddyCurrentFoam'

# Init dictionaries as we are using meshes
data[set] = dict()

R[set] = dict()
Z[set] = dict()

F[set] = dict()

E[set] = dict()
D[set] = dict()
N[set] = dict()

meshes = ['0.125', '0.250', '0.375', '0.500', '0.750', '1.000', '1.500', '2.000']

# Read data
for mesh in meshes:

    data[set][mesh] = np.genfromtxt(fileGetPath(baseName+set+'_'+mesh+'.dat'), comments='#')

    R[set][mesh]    = data[set][mesh][:,0].reshape(nr,nz)
    Z[set][mesh]    = data[set][mesh][:,2].reshape(nr,nz)

    F[set][mesh]    = [np.zeros(R[set][mesh].shape) for i in range(3)]
    F[set][mesh][0] = data[set][mesh][:,3].reshape(nr,nz)
    F[set][mesh][1] = data[set][mesh][:,4].reshape(nr,nz)
    F[set][mesh][2] = data[set][mesh][:,5].reshape(nr,nz)

    # Scale to mm
    R[set][mesh]    = 1000.0 * (R[set][mesh])
    Z[set][mesh]    = 1000.0 * (Z[set][mesh] - 0.03)

    # Flip y-axis
    F[set][mesh][1] = -F[set][mesh][1]

# Calculate errors with repect to last mesh in list
for mesh in meshes[:-1]:

    E[set][mesh]    = [np.zeros(R[set][mesh].shape) for i in range(3)]

    for i in range(3):

        E[set][mesh][i] = abs(F[set][meshes[-1]][i] - F[set][mesh][i])

# Init norms
for key in norms.keys():

    N[set][key] = dict()

    for i in range(3):

        N[set][key][i] = dict()

# Calculate norms and global maximum
Dmax = 0.0
Nmax = [ 0.0 for i in range(3)]
for mesh in meshes[:-1]:

    D[set][mesh] = 1.0/float(mesh)
    Dmax = max(Dmax, D[set][mesh])

    for key in norms.keys():

        for i in range(3):

            N[set][key][i][mesh] = norms[key](F[set][meshes[-1]][i], E[set][mesh][i])
            Nmax[i] = max(Nmax[i], N[set][key][i][mesh])

# Scale norms and mesh delta
for mesh in meshes[:-1]:

    D[set][mesh] /= Dmax

    for key in norms.keys():

        for i in range(3):

            N[set][key][i][mesh] /= Nmax[i]

#print [ i for k, i in sorted(N[set]['inf'][0].iteritems())]
#print [ i for k, i in sorted(N[set]['1'][0].iteritems())]
#print [ i for k, i in sorted(N[set]['2'][0].iteritems())]
#print [ i for k, i in sorted(D[set].iteritems())]

# --------------------------------------------------------------------------- #
# --- Plot settings --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

latexify(fontsize=fontsize, fontfamily=fontfamily)

labelR = '$r ~ [\mathrm{mm}]$'
labelZ = '$z ~ [\mathrm{mm}]$'

labelEinf = '$\mathrm{log}\|E\|_{\infty}$'
labelE1 = '$\mathrm{log}\|E\|_{1}$'
labelE2 = '$\mathrm{log}\|E\|_{2}$'

levelsFr = np.linspace(-0.30,-0.05,6)
levelsFa = np.linspace( 0.20, 2.00,10)
levelsFz = np.linspace(-0.40,-0.05,8)
levelsF  = levelsFa

colors = 'black'

plots = dict()

# --------------------------------------------------------------------------- #
# --- Error plots ----------------------------------------------------------- #
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

        axe.set_xlim([1,5e-2])
        axe.set_ylim([1e-4,1])

        axe.set_xscale('log')
        axe.set_yscale('log')

        def ele(f, a, elements, name):

            set = 'EddyCurrentFoam'

            elements[name] = a.plot([ i for k, i in sorted(D[set].iteritems())], [ i for k, i in sorted(N[set]['inf'][0].iteritems())], label=labelEinf)
            elements[name] = a.plot([ i for k, i in sorted(D[set].iteritems())], [ i for k, i in sorted(N[set]['1'][0].iteritems())], label=labelE1)
            elements[name] = a.plot([ i for k, i in sorted(D[set].iteritems())], [ i for k, i in sorted(N[set]['2'][0].iteritems())], label=labelE2)

            a.legend()

        ele(fig, axe, elem, 'norms')

    axe(fig, axes, 'Fr')

    fig.savefig(fileGetPath(baseName+name+'.pdf'), bbox_inches="tight")

fig(plots, 'ErrorFr')



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

        axe.set_xlim([1,5e-2])
        axe.set_ylim([1e-4,1])

        axe.set_xscale('log')
        axe.set_yscale('log')

        def ele(f, a, elements, name):

            set = 'EddyCurrentFoam'

            elements[name] = a.plot([ i for k, i in sorted(D[set].iteritems())], [ i for k, i in sorted(N[set]['inf'][1].iteritems())], label=labelEinf)
            elements[name] = a.plot([ i for k, i in sorted(D[set].iteritems())], [ i for k, i in sorted(N[set]['1'][1].iteritems())], label=labelE1)
            elements[name] = a.plot([ i for k, i in sorted(D[set].iteritems())], [ i for k, i in sorted(N[set]['2'][1].iteritems())], label=labelE2)

            a.legend()

        ele(fig, axe, elem, 'norms')

    axe(fig, axes, 'Fa')

    fig.savefig(fileGetPath(baseName+name+'.pdf'), bbox_inches="tight")

fig(plots, 'ErrorFa')



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

        axe.set_xlim([1,5e-2])
        axe.set_ylim([1e-4,1])

        axe.set_xscale('log')
        axe.set_yscale('log')

        def ele(f, a, elements, name):

            set = 'EddyCurrentFoam'

            elements[name] = a.plot([ i for k, i in sorted(D[set].iteritems())], [ i for k, i in sorted(N[set]['inf'][2].iteritems())], label=labelEinf)
            elements[name] = a.plot([ i for k, i in sorted(D[set].iteritems())], [ i for k, i in sorted(N[set]['1'][2].iteritems())], label=labelE1)
            elements[name] = a.plot([ i for k, i in sorted(D[set].iteritems())], [ i for k, i in sorted(N[set]['2'][2].iteritems())], label=labelE2)

            a.legend()

        ele(fig, axe, elem, 'norms')

    axe(fig, axes, 'Fz')

    fig.savefig(fileGetPath(baseName+name+'.pdf'), bbox_inches="tight")

fig(plots, 'ErrorFz')

# --------------------------------------------------------------------------- #
# --- Contour plots --------------------------------------------------------- #
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

            set = 'Analytical'

            elements[name] = a.contour(
                R[set], Z[set], F[set][1],
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

            set = 'Opera3D'

            elements[name] = a.contour(
                R[set], Z[set], F[set][0],
                levels=levelsFr, colors=colors, linestyles='dashed')

        ele(fig, axe, elem, 'cOpera3D')

        def ele(f, a, elements, name):

            set = 'EddyCurrentFoam'
            mesh = '1.000'

            elements[name] = a.contour(
                R[set][mesh], Z[set][mesh], F[set][mesh][0],
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

            set = 'Analytical'

            elements[name] = a.contour(
                R[set], Z[set], F[set][1],
                levels=levelsFa, colors=colors, linestyles='dotted')

        ele(fig, axe, elem, 'cAnalytical')

        def ele(f, a, elements, name):

            set = 'Opera3D'

            elements[name] = a.contour(
                R[set], Z[set], F[set][1],
                levels=levelsFr, colors=colors, linestyles='dashed')

        ele(fig, axe, elem, 'cOpera3D')

        def ele(f, a, elements, name):

            set = 'EddyCurrentFoam'
            mesh = '1.000'

            elements[name] = a.contour(
                R[set][mesh], Z[set][mesh], F[set][mesh][1],
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

            set = 'Opera3D'

            elements[name] = a.contour(
                R[set], Z[set], F[set][2],
                levels=levelsFr, colors=colors, linestyles='dashed')

        ele(fig, axe, elem, 'cOpera3D')

        def ele(f, a, elements, name):

            set = 'EddyCurrentFoam'
            mesh = '1.000'

            elements[name] = a.contour(
                R[set][mesh], Z[set][mesh], F[set][mesh][2],
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

            set = 'Opera3D'

            magF = (F[set][0]**2 + F[set][1]**2 + F[set][2]**2)**0.5

            elements[name] = a.contour(
                R[set], Z[set], magF,
                levels=levelsF, colors=colors, linestyles='dashed')

        ele(fig, axe, elem, 'cOpera3D')

        def ele(f, a, elements, name):

            set = 'EddyCurrentFoam'
            mesh = '1.000'

            magF = (F[set][mesh][0]**2 + F[set][mesh][1]**2 + F[set][mesh][2]**2)**0.5

            elements[name] = a.contour(
                R[set][mesh], Z[set][mesh], magF,
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
