#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
# Pascal Beckstein (p.beckstein@hzdr.de)

from __future__ import nested_scopes, generators, division, absolute_import
from __future__ import with_statement, print_function, unicode_literals

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

__name__
__path__ = os.path.realpath(__file__)
__base__ = os.path.basename(__path__)
__dir__ = os.path.dirname(__path__)
__head__ = os.path.splitext(__base__)[0]

sys.path.append(os.environ["FOAM_USER_TOOLS"] + "/" + "python")

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

fontsize   = 16
fontfamily = "serif"

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
names = dict()

# --------------------------------------------------------------------------- #

def readdata(set):

    cases = ["ortho", "nonortho"]
    meshes = ["1.000"]
    lines = ["x1", "y1", "y2", "z1"]
    frequencies = ["1000", "10000", "100000"]

    setData = dict()
    setNames = dict()

    for case in cases:

        setData[case] = dict()
        setNames[case] = dict()

        for mesh in meshes:

            setData[case][mesh] = dict()
            setNames[case][mesh] = dict()

            for freq in frequencies:

                setData[case][mesh][freq] = dict()
                setNames[case][mesh][freq] = dict()

                for line in lines:

                    fileName = __dir__ + "/data_" + set \
                             + "_" + case + "_m" + mesh \
                             + "_f" + freq + "_line_" + line + ".dat"

                    if os.path.isfile(fileName):

                        print("Reading file: " + fileName)

                        setData[case][mesh][freq][line] = \
                            np.genfromtxt(fileName, comments='#', names=True)
                        setNames[case][mesh][freq][line] = \
                            setData[case][mesh][freq][line].dtype.names

    return setData, setNames

# --------------------------------------------------------------------------- #

set = "Opera3D"

data[set], names[set] = readdata(set)

# --------------------------------------------------------------------------- #

set = "eddyCurrentFoam"

data[set], names[set] = readdata(set)

# --------------------------------------------------------------------------- #
# --- Plot settings --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

latex.latexify(fontsize=fontsize, fontfamily=fontfamily, locale="de_DE.utf8")

hzdr.colors()

plots = dict()

# --------------------------------------------------------------------------- #
# --- Test ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

def fig(p, name):

    p[name] = {"fig": plt.figure(), "axs": dict()}
    f = p[name]

    fig = f["fig"]
    axs = f["axs"]

    def ax(f, axs, name):

        axs[name] = fig.add_subplot(111)
        ax = axs[name]

        ax.minorticks_on()

        ax.set_xlim([-75,75])
        #ax.set_ylim([0,30])

        #ax.set_xlabel(labelAxisR)
        #ax.set_ylabel(labelAxisZ)

        #ax.set_aspect("equal")

        opData = data["Opera3D"]["ortho"]["1.000"]["1000"]["y2"]
        ofData = data["eddyCurrentFoam"]["ortho"]["1.000"]["1000"]["y2"]

        complexNames = ["Re", "Im"]
        dirNames = ["x", "y", "z"]

        colors = dict()

        colors["Re"] = ["hzdr-orange", "hzdr-yellow", "hzdr-green"]
        colors["Im"] = ["hzdr-blue", "hzdr-purple", "hzdr-red"]

        markers = dict()

        markers["Re"] = ["o", "v", "^"]
        markers["Im"] = ["s", "d", "*"]

        # j in [10^6 A/m^-2]
        for c in complexNames:

            for i, d in enumerate(dirNames):

                var = "j" + c + "_" + d
                color = colors[c][i]
                marker = markers[c][i]
                label = r"${\boldsymbol{j}_" + d + r"}_{\,\scriptstyle\mathfrak{" + c + r"}}$"

                ax.plot(opData["y"], opData[var], color=color, linestyle="--")

                ax.plot(1e+3*ofData["y"], 1e-6*ofData[var], color=color, linestyle="-",
                        marker=marker, markevery=5, markersize=5,
                        markeredgecolor=color, markerfacecolor=color,
                        label=label)

        ax.legend(bbox_to_anchor=(0.0, 1.1, 1.0, .1), loc="upper center", ncol=3, mode="expand", borderaxespad=0.)

        #[10^-2 T]
        #ax.plot(opData["y"], 1e+2*opData["BRe_y"], label="Opera3D")
        #ax.plot(1e+3*ofData["y"], 1e+2*ofData["BRe_y"], label="eddyCurrentFoam")

        ##[10^4 N/m^3]
        #ax.plot(opData["y"], 1e+2*opData["F_y"], label="Opera3D")
        #ax.plot(1e+3*ofData["y"], 1e-4*ofData["F_y"], label="eddyCurrentFoam")

        ax.legend(bbox_to_anchor=(0.0, 1.1, 1.0, .1), loc="upper center", ncol=3, mode="expand", borderaxespad=0.)

        #set = "Analytical"
        #c = ax.contour(R[set], Z[set], F[set][1],
                       #levels=levels[1])
        #cl = ax.clabel(c, c.levels[0::2],
                       #inline=True, fmt="%g", fontsize=fontsize)

        #[l.set_bbox(dict(facecolor="white", edgecolor="none", pad=2)) for l in cl]
        #[l.set_text('{:n}'.format(float(l.get_text()))) for l in cl]

    ax(fig, axs, "test")

    #fig.set_size_inches(sizeCompX, sizeCompY)
    #fig.savefig(__dir__+"/"+baseName+name+".pdf", bbox_inches="tight")
    fig.savefig(__dir__+ "/" + name + ".pdf", bbox_inches="tight")

fig(plots, "test")

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
