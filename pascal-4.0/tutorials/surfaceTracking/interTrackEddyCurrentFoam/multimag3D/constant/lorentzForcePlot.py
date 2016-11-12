#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

csb = os.path.basename(os.path.realpath(sys.argv[0]))
csd = os.path.dirname(os.path.realpath(sys.argv[0]))
csn = os.path.splitext(csb)[0]

sys.path.append(os.environ["FOAM_USER_TOOLS"]+"/python")

import math as m
import numpy as np

import matplotlib.pyplot as plt

import plotTools.latex as latex
import plotTools.hzdr as hzdr

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

nr = 61
nz = 61

fontsize   = 14
fontfamily = "serif"

print "nr, nz     : ", nr, ",", nz
print "fontsize   : ", fontsize
print "fontfamily : ", fontfamily

baseName = "lorentzForce"

# --------------------------------------------------------------------------- #
# --- Functions ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def norm_inf(r, e): return np.max(np.absolute(e)) / np.max(np.absolute(r))
def norm_1(r, e): return np.sum(np.absolute(e)) / np.sum(np.absolute(r))
def norm_2(r, e): return (np.sum(e**2) / np.sum(r**2))**0.5

norms = {"inf": norm_inf, "1": norm_1, "2": norm_2}

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

set = "Analytical"

# Read data
if True:

    data[set] = np.genfromtxt(csd+"/"+baseName+set+".dat", comments="#")

    R[set]    = data[set][:,0].reshape(nr,nz)
    Z[set]    = data[set][:,1].reshape(nr,nz)

    F[set]    = [np.zeros(R[set].shape) for i in range(3)]
    F[set][1] = data[set][:,2].reshape(nr,nz)

    # Scale to mm
    R[set] = 1000.0 * R[set]
    Z[set] = 1000.0 * Z[set]

# --------------------------------------------------------------------------- #

set = "Opera3D"

# Init dictionaries as we are using meshes
data[set] = dict()

R[set] = dict()
Z[set] = dict()

F[set] = dict()

E[set] = dict()
D[set] = dict()
N[set] = dict()

meshes = ["coarse", "fine"]

# Read data
for mesh in meshes:

    data[set][mesh] = np.genfromtxt(csd+"/"+baseName+set+"_"+mesh+".dat", comments="#")

    R[set][mesh]    = data[set][mesh][:,0].reshape(nr,nz)
    Z[set][mesh]    = data[set][mesh][:,1].reshape(nr,nz)

    F[set][mesh]    = [np.zeros(R[set][mesh].shape) for i in range(3)]
    F[set][mesh][0] = data[set][mesh][:,2].reshape(nr,nz)
    F[set][mesh][1] = data[set][mesh][:,3].reshape(nr,nz)
    F[set][mesh][2] = data[set][mesh][:,4].reshape(nr,nz)

# --------------------------------------------------------------------------- #

set = "EddyCurrentFoam"

# Init dictionaries as we are using meshes
data[set] = dict()

R[set] = dict()
Z[set] = dict()

F[set] = dict()

E[set] = dict()
D[set] = dict()
N[set] = dict()

meshes = ["0.125", "0.250", "0.375", "0.500", "0.750", "1.000", "1.500", "2.000"]

# Read data
for mesh in meshes:

    data[set][mesh] = np.genfromtxt(csd+"/"+baseName+set+"_"+mesh+".dat", comments="#")

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
for norm in norms.keys():

    N[set][norm] = dict()

    for i in range(3):

        N[set][norm][i] = dict()

# Calculate norms and global maximum
Dmax = 0.0
Nmax = [ 0.0 for i in range(3)]
for mesh in meshes[:-1]:

    D[set][mesh] = 1.0/float(mesh)
    Dmax = max(Dmax, D[set][mesh])

    for norm in norms.keys():

        for i in range(3):

            N[set][norm][i][mesh] = norms[norm](F[set][meshes[-1]][i], E[set][mesh][i])
            Nmax[i] = max(Nmax[i], N[set][norm][i][mesh])

# Scale norms and mesh delta
for mesh in meshes[:-1]:

    D[set][mesh] /= Dmax

    for norm in norms.keys():

        for i in range(3):

            N[set][norm][i][mesh] /= Nmax[i]

# --------------------------------------------------------------------------- #

set = "EddyCurrentFoam_space"

# Init dictionaries as we are using meshes
data[set] = dict()

R[set] = dict()
Z[set] = dict()

F[set] = dict()

E[set] = dict()
D[set] = dict()
N[set] = dict()

meshes = ["00.300", "00.375", "00.500", "00.750", "01.000", "02.000", "04.000", "08.000", "16.000", "32.000"]

# Read data
for mesh in meshes:

    data[set][mesh] = np.genfromtxt(csd+"/"+baseName+set+"_"+mesh+".dat", comments="#")

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
for norm in norms.keys():

    N[set][norm] = dict()

    for i in range(3):

        N[set][norm][i] = dict()

# Calculate norms and global maximum
#Dmax = 0.0
Nmax = [ 0.0 for i in range(3)]
for mesh in meshes[:-1]:

    D[set][mesh] = 1.0/float(mesh)
    D[set][mesh] = float(mesh)
    Dmax = max(Dmax, D[set][mesh])

    for norm in norms.keys():

        for i in range(3):

            N[set][norm][i][mesh] = norms[norm](F[set][meshes[-1]][i], E[set][mesh][i])
            Nmax[i] = max(Nmax[i], N[set][norm][i][mesh])

# Scale norms and space size
for mesh in meshes[:-1]:

    #D[set][mesh] /= Dmax
    D[set][mesh] *= 4.0

    for norm in norms.keys():

        for i in range(3):

            N[set][norm][i][mesh] /= Nmax[i]

# --------------------------------------------------------------------------- #
# --- Plot settings --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

latex.latexify(fontsize=fontsize, fontfamily=fontfamily)

hzdr.colors()

labelAxisR = r"$r ~ [\mathrm{mm}]$"
labelAxisZ = r"$z ~ [\mathrm{mm}]$"

labelAxisE = r"$\mathrm{log}\|E\|$"
labelAxisD = r"$\mathrm{log}(\triangle x /\triangle x_{\mathrm{ref}})$"
labelAxisS = r"$\mathrm{log}(R_\mathrm{0}/R_\mathrm{\bf j})$"

labelE = dict()
labelE["O1"] = r"$\mathcal{O}(\triangle x /\triangle x_{\mathrm{ref}})$"
labelE["O2"] = r"$\mathcal{O}((\triangle x /\triangle x_{\mathrm{ref}})^2)$"
labelE["Os3"] = r"$\mathcal{O}((R_\mathrm{0}/R_\mathrm{\bf j})^{-3})$"
labelE["OsConst"] = r"$R_\mathrm{0}/R_\mathrm{\bf j} = r$"
labelE["inf"] = r"$\|E\|_{\infty}$"
labelE["1"] = r"$\|E\|_{1}$"
labelE["2"] = r"$\|E\|_{2}$"

markerE = dict()
markerE["O1"] = ""
markerE["O2"] = ""
markerE["Os3"] = ""
markerE["inf"] = "."
markerE["1"] = "x"
markerE["2"] = "+"

levels = [np.linspace(-0.30,-0.05,6),
          np.linspace( 0.20, 2.00,10),
          np.linspace(-0.40,-0.05,8)]

#colors = "black"

plots = dict()

# --------------------------------------------------------------------------- #
# --- Error plots ----------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([1,5e-2])
        ax.set_ylim([1e-4,1])

        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.set_xlabel(labelAxisD)
        ax.set_ylabel(labelAxisE)

        d = np.linspace(1,1e-2,100)
        ax.plot(d, 1.8*d, label=labelE["O1"],
                marker=markerE["O1"], linestyle="dotted", color="black")
        ax.plot(d, 0.2*d**2.0, label=labelE["O2"],
                marker=markerE["O2"], linestyle="dashed", color="black")

        set = "EddyCurrentFoam"

        d = np.array([ i for k, i in sorted(D[set].iteritems())])
        for norm in sorted(norms.keys()):
            n = np.array([ i for k, i in sorted(N[set][norm][0].iteritems())])
            ax.plot(d, n, label=labelE[norm], marker=markerE[norm])

        ax.legend(loc="lower left")

    ax(fig, axs, "Fr")

    fig.savefig(csd+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "ErrorFr")



def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([1,5e-2])
        ax.set_ylim([1e-4,1])

        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.set_xlabel(labelAxisD)
        ax.set_ylabel(labelAxisE)

        d = np.linspace(1,1e-2,100)
        ax.plot(d, 1.8*d, label=labelE["O1"],
                marker=markerE["O1"], linestyle="dotted", color="black")
        ax.plot(d, 0.2*d**2.0, label=labelE["O2"],
                marker=markerE["O2"], linestyle="dashed", color="black")

        set = "EddyCurrentFoam"

        d = np.array([ i for k, i in sorted(D[set].iteritems())])
        for norm in sorted(norms.keys()):
            n = np.array([ i for k, i in sorted(N[set][norm][1].iteritems())])
            ax.plot(d, n, label=labelE[norm], marker=markerE[norm])

        ax.legend(loc="lower left")

    ax(fig, axs, "Fa")

    fig.savefig(csd+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "ErrorFa")



def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([1,5e-2])
        ax.set_ylim([1e-4,1])

        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.set_xlabel(labelAxisD)
        ax.set_ylabel(labelAxisE)

        d = np.linspace(1,1e-2,100)
        ax.plot(d, 1.8*d, label=labelE["O1"],
                marker=markerE["O1"], linestyle="dotted", color="black")
        ax.plot(d, 0.2*d**2.0, label=labelE["O2"],
                marker=markerE["O2"], linestyle="dashed", color="black")

        set = "EddyCurrentFoam"

        d = np.array([ i for k, i in sorted(D[set].iteritems())])
        for norm in sorted(norms.keys()):
            n = np.array([ i for k, i in sorted(N[set][norm][2].iteritems())])
            ax.plot(d, n, label=labelE[norm], marker=markerE[norm])

        ax.legend(loc="lower left")

    ax(fig, axs, "Fz")

    fig.savefig(csd+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "ErrorFz")



def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([1,5e-2])
        ax.set_ylim([1e-4,1])

        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.set_xlabel(labelAxisD)
        ax.set_ylabel(labelAxisE)

        d = np.linspace(1,1e-2,100)
        ax.plot(d, 1.8*d, label=labelE["O1"],
                marker=markerE["O1"], linestyle="dotted", color="black")
        ax.plot(d, 0.2*d**2.0, label=labelE["O2"],
                marker=markerE["O2"], linestyle="dashed", color="black")

        set = "EddyCurrentFoam"

        d = np.array([ i for k, i in sorted(D[set].iteritems())])
        magnMax = 0.0
        magn = dict()
        for norm in sorted(norms.keys()):
            n0 = np.array([ i for k, i in sorted(N[set][norm][0].iteritems())])
            n1 = np.array([ i for k, i in sorted(N[set][norm][1].iteritems())])
            n2 = np.array([ i for k, i in sorted(N[set][norm][2].iteritems())])
            magn[norm] = (n0**2 + n1**2 + n2**2)**0.5
            for m in magn[norm]:
                magnMax = max(magnMax, m)

        for norm in sorted(norms.keys()):
            n = magn[norm]/magnMax
            ax.plot(d, n, label=labelE[norm], marker=markerE[norm])

        ax.legend(loc="lower left")

    ax(fig, axs, "F")

    fig.savefig(csd+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "ErrorF")

# --------------------------------------------------------------------------- #
# --- Space error plots ----------------------------------------------------- #
# --------------------------------------------------------------------------- #

def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([1,1e2])
        ax.set_ylim([1e-5,1])

        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.set_xlabel(labelAxisS)
        ax.set_ylabel(labelAxisE)

        for n in [1, 2, 3, 4]:

            c = 10.0**(n/3.0)

            if n==1: y = 5e-4
            if n==2: y = 5e-5
            if n==3: y = 5e-2
            if n==4: y = 5e-3

            if n==1:
                ax.plot((c, c), (1e-5, 1), label=labelE["OsConst"],
                        linestyle="dotted", color="black")
            else:
                ax.plot((c, c), (1e-5, 1),
                        linestyle="dotted", color="black")

            ax.text(c, y, "$r = {:.1f}$".format(c),
                    horizontalalignment="center", verticalalignment="center",
                    bbox=dict(facecolor="white", edgecolor="none", pad=2))

        d = np.linspace(1,100,100)
        ax.plot(d, 1.0*d**-3.0, label=labelE["Os3"],
                marker=markerE["Os3"], linestyle="dashed", color="black")

        set = "EddyCurrentFoam_space"

        d = np.array([ i for k, i in sorted(D[set].iteritems())])
        for norm in sorted(norms.keys()):
            n = np.array([ i for k, i in sorted(N[set][norm][0].iteritems())])
            ax.plot(d, n, label=labelE[norm], marker=markerE[norm])

        ax.legend(loc="upper right")

    ax(fig, axs, "Fr")

    fig.savefig(csd+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "SpaceErrorFr")



def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([1,1e2])
        ax.set_ylim([1e-5,1])

        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.set_xlabel(labelAxisS)
        ax.set_ylabel(labelAxisE)

        for n in [1, 2, 3, 4]:

            c = 10.0**(n/3.0)

            if n==1: y = 5e-4
            if n==2: y = 5e-5
            if n==3: y = 5e-2
            if n==4: y = 5e-3

            if n==1:
                ax.plot((c, c), (1e-5, 1), label=labelE["OsConst"],
                        linestyle="dotted", color="black")
            else:
                ax.plot((c, c), (1e-5, 1),
                        linestyle="dotted", color="black")

            ax.text(c, y, "$r = {:.1f}$".format(c),
                    horizontalalignment="center", verticalalignment="center",
                    bbox=dict(facecolor="white", edgecolor="none", pad=2))

        d = np.linspace(1,100,100)
        ax.plot(d, 1.0*d**-3.0, label=labelE["Os3"],
                marker=markerE["Os3"], linestyle="dashed", color="black")

        set = "EddyCurrentFoam_space"

        d = np.array([ i for k, i in sorted(D[set].iteritems())])
        for norm in sorted(norms.keys()):
            n = np.array([ i for k, i in sorted(N[set][norm][1].iteritems())])
            ax.plot(d, n, label=labelE[norm], marker=markerE[norm])

        ax.legend(loc="upper right")

    ax(fig, axs, "Fa")

    fig.savefig(csd+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "SpaceErrorFa")



def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([1,1e2])
        ax.set_ylim([1e-5,1])

        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.set_xlabel(labelAxisS)
        ax.set_ylabel(labelAxisE)

        for n in [1, 2, 3, 4]:

            c = 10.0**(n/3.0)

            if n==1: y = 5e-4
            if n==2: y = 5e-5
            if n==3: y = 5e-2
            if n==4: y = 5e-3

            if n==1:
                ax.plot((c, c), (1e-5, 1), label=labelE["OsConst"],
                        linestyle="dotted", color="black")
            else:
                ax.plot((c, c), (1e-5, 1),
                        linestyle="dotted", color="black")

            ax.text(c, y, "$r = {:.1f}$".format(c),
                    horizontalalignment="center", verticalalignment="center",
                    bbox=dict(facecolor="white", edgecolor="none", pad=2))

        d = np.linspace(1,100,100)
        ax.plot(d, 1.0*d**-3.0, label=labelE["Os3"],
                marker=markerE["Os3"], linestyle="dashed", color="black")

        set = "EddyCurrentFoam_space"

        d = np.array([ i for k, i in sorted(D[set].iteritems())])
        for norm in sorted(norms.keys()):
            n = np.array([ i for k, i in sorted(N[set][norm][2].iteritems())])
            ax.plot(d, n, label=labelE[norm], marker=markerE[norm])

        ax.legend(loc="upper right")

    ax(fig, axs, "Fz")

    fig.savefig(csd+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "SpaceErrorFz")



def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([1,1e2])
        ax.set_ylim([1e-5,1])

        ax.set_xscale("log")
        ax.set_yscale("log")

        ax.set_xlabel(labelAxisS)
        ax.set_ylabel(labelAxisE)

        for n in [1, 2, 3, 4]:

            c = 10.0**(n/3.0)

            if n==1: y = 5e-4
            if n==2: y = 5e-5
            if n==3: y = 5e-2
            if n==4: y = 5e-3

            if n==1:
                ax.plot((c, c), (1e-5, 1), label=labelE["OsConst"],
                        linestyle="dotted", color="black")
            else:
                ax.plot((c, c), (1e-5, 1),
                        linestyle="dotted", color="black")

            ax.text(c, y, "$r = {:.1f}$".format(c),
                    horizontalalignment="center", verticalalignment="center",
                    bbox=dict(facecolor="white", edgecolor="none", pad=2))

        d = np.linspace(1,100,100)
        ax.plot(d, 1.0*d**-3.0, label=labelE["Os3"],
                marker=markerE["Os3"], linestyle="dashed", color="black")

        set = "EddyCurrentFoam_space"

        d = np.array([ i for k, i in sorted(D[set].iteritems())])
        magnMax = 0.0
        magn = dict()
        for norm in sorted(norms.keys()):
            n0 = np.array([ i for k, i in sorted(N[set][norm][0].iteritems())])
            n1 = np.array([ i for k, i in sorted(N[set][norm][1].iteritems())])
            n2 = np.array([ i for k, i in sorted(N[set][norm][2].iteritems())])
            magn[norm] = (n0**2 + n1**2 + n2**2)**0.5
            for m in magn[norm]:
                magnMax = max(magnMax, m)

        for norm in sorted(norms.keys()):
            n = magn[norm]/magnMax
            ax.plot(d, n, label=labelE[norm], marker=markerE[norm])

        ax.legend(loc="upper right")

    ax(fig, axs, "F")

    fig.savefig(csd+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "SpaceErrorF")

# --------------------------------------------------------------------------- #
# --- Contour plots --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([0,30])
        ax.set_ylim([0,30])

        ax.set_xlabel(labelAxisR)
        ax.set_ylabel(labelAxisZ)

        ax.set_aspect("equal")

        set = "Analytical"
        c = ax.contour(R[set], Z[set], F[set][1],
                       levels=levels[1])
        cl = ax.clabel(c, c.levels[0::2],
                       inline=True, fmt="%g", fontsize=fontsize)

    ax(fig, axs, "Fa")

    fig.savefig(csd+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "Analytical")



def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([0,30])
        ax.set_ylim([0,30])

        ax.set_xlabel(labelAxisR)
        ax.set_ylabel(labelAxisZ)

        ax.set_aspect("equal")

        set = "Opera3D"
        mesh = "fine"
        c = ax.contour(R[set][mesh], Z[set][mesh], F[set][mesh][0],
                       levels=levels[0], linestyles="dashed")

        set = "EddyCurrentFoam"
        mesh = "1.000"
        c = ax.contour(R[set][mesh], Z[set][mesh], F[set][mesh][0],
                       levels=levels[0], linestyles="solid")
        cl = ax.clabel(c, c.levels[0::2],
                       inline=True, fmt="%g", fontsize=fontsize)

        [l.set_bbox(dict(facecolor="white", edgecolor="none", pad=2)) for l in cl]

    ax(fig, axs,"Fr")

    fig.savefig(csd+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "ComparisonFr")



def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([0,30])
        ax.set_ylim([0,30])

        ax.set_xlabel(labelAxisR)
        ax.set_ylabel(labelAxisZ)

        ax.set_aspect("equal")

        set = "Analytical"
        ax.contour(R[set], Z[set], F[set][1],
                   levels=levels[1], linestyles="dotted")

        set = "Opera3D"
        mesh = "fine"
        c = ax.contour(R[set][mesh], Z[set][mesh], F[set][mesh][1],
                       levels=levels[1], linestyles="dashed")

        set = "EddyCurrentFoam"
        mesh = "1.000"
        c = ax.contour(R[set][mesh], Z[set][mesh], F[set][mesh][1],
                       levels=levels[1], linestyles="solid")
        cl = ax.clabel(c, c.levels[0::2],
                       inline=True, fmt="%g", fontsize=fontsize)

        [l.set_bbox(dict(facecolor="white", edgecolor="none", pad=2)) for l in cl]

    ax(fig, axs, "Fa")

    fig.savefig(csd+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "ComparisonFa")



def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([0,30])
        ax.set_ylim([0,30])

        ax.set_xlabel(labelAxisR)
        ax.set_ylabel(labelAxisZ)

        ax.set_aspect("equal")

        set = "Opera3D"
        mesh = "fine"
        c = ax.contour(R[set][mesh], Z[set][mesh], F[set][mesh][2],
                       levels=levels[2], linestyles="dashed")

        set = "EddyCurrentFoam"
        mesh = "1.000"
        c = ax.contour(R[set][mesh], Z[set][mesh], F[set][mesh][2],
                       levels=levels[2], linestyles="solid")
        cl = ax.clabel(c, c.levels[0::2],
                       inline=True, fmt="%g", fontsize=fontsize)

        [l.set_bbox(dict(facecolor="white", edgecolor="none", pad=2)) for l in cl]

    ax(fig, axs, "Fz")

    fig.savefig(csd+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "ComparisonFz")



def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.set_xlim([0,30])
        ax.set_ylim([0,30])

        ax.set_xlabel(labelAxisR)
        ax.set_ylabel(labelAxisZ)

        ax.set_aspect("equal")

        set = "Opera3D"
        mesh = "fine"
        magF = (F[set][mesh][0]**2 + F[set][mesh][1]**2 + F[set][mesh][2]**2)**0.5
        c = ax.contour(R[set][mesh], Z[set][mesh], magF,
                       levels=levels[1], linestyles="dashed")

        set = "EddyCurrentFoam"
        mesh = "1.000"
        magF = (F[set][mesh][0]**2 + F[set][mesh][1]**2 + F[set][mesh][2]**2)**0.5
        c = ax.contour(R[set][mesh], Z[set][mesh], magF,
                       levels=levels[1], linestyles="solid")
        cl = ax.clabel(c, c.levels[0::2],
                       inline=True, fmt="%g", fontsize=fontsize)

        [l.set_bbox(dict(facecolor="white", edgecolor="none", pad=2)) for l in cl]

    ax(fig, axs, "F")

    fig.savefig(csd+"/"+baseName+name+".pdf", bbox_inches="tight")

fig(plots, "ComparisonF")

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
