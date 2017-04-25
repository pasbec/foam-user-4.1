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
scales = dict()

# --------------------------------------------------------------------------- #

def readdata(set):

    cases = ["ortho", "nonortho"]
    meshes = ["0.125", "0.250", "0.375", "0.500", "0.750", "1.000", "1.500", "2.000", "2.500", "3.000"]
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

scales[set] = dict()
scales[set]["x"] = 1e-3
scales[set]["y"] = 1e-3
scales[set]["z"] = 1e-3
scales[set]["jRe_x"] = 1e+6
scales[set]["jRe_y"] = 1e+6
scales[set]["jRe_z"] = 1e+6
scales[set]["jIm_x"] = 1e+6
scales[set]["jIm_y"] = 1e+6
scales[set]["jIm_z"] = 1e+6
scales[set]["BRe_x"] = 1.0
scales[set]["BRe_y"] = 1.0
scales[set]["BRe_z"] = 1.0
scales[set]["BIm_x"] = 1.0
scales[set]["BIm_y"] = 1.0
scales[set]["BIm_z"] = 1.0
scales[set]["F_x"] = 1e+6
scales[set]["F_y"] = 1e+6
scales[set]["F_z"] = 1e+6

# --------------------------------------------------------------------------- #

set = "eddyCurrentFoam"

data[set], names[set] = readdata(set)

scales[set] = dict()
scales[set]["x"] = 1.0
scales[set]["y"] = 1.0
scales[set]["z"] = 1.0
scales[set]["jRe_x"] = 1.0
scales[set]["jRe_y"] = 1.0
scales[set]["jRe_z"] = 1.0
scales[set]["jIm_x"] = 1.0
scales[set]["jIm_y"] = 1.0
scales[set]["jIm_z"] = 1.0
scales[set]["BRe_x"] = 1.0
scales[set]["BRe_y"] = 1.0
scales[set]["BRe_z"] = 1.0
scales[set]["BIm_x"] = 1.0
scales[set]["BIm_y"] = 1.0
scales[set]["BIm_z"] = 1.0
scales[set]["F_x"] = 1.0
scales[set]["F_y"] = 1.0
scales[set]["F_z"] = 1.0
# VReGrad_x VReGrad_y VReGrad_z VImGrad_x VImGrad_y VImGrad_z VRe VIm sigma mur

# --------------------------------------------------------------------------- #
# --- Plot settings --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

latex.latexify(fontsize=fontsize, fontfamily=fontfamily, locale="de_DE.utf8")

hzdr.colors()

alabels = dict()
alabels["x"] = r"$x ~ [\mathrm{mm}]$"
alabels["y"] = r"$y ~ [\mathrm{mm}]$"
alabels["z"] = r"$z ~ [\mathrm{mm}]$"

labels = dict()
labels["jRe_x"] = r"${j_x}_{\,\scriptstyle\mathfrak{Re}}$"
labels["jRe_y"] = r"${j_y}_{\,\scriptstyle\mathfrak{Re}}$"
labels["jRe_z"] = r"${j_z}_{\,\scriptstyle\mathfrak{Re}}$"
labels["jIm_x"] = r"${j_x}_{\,\scriptstyle\mathfrak{Im}}$"
labels["jIm_y"] = r"${j_y}_{\,\scriptstyle\mathfrak{Im}}$"
labels["jIm_z"] = r"${j_z}_{\,\scriptstyle\mathfrak{Im}}$"
labels["BRe_x"] = r"${B_x}_{\,\scriptstyle\mathfrak{Re}}$"
labels["BRe_y"] = r"${B_y}_{\,\scriptstyle\mathfrak{Re}}$"
labels["BRe_z"] = r"${B_z}_{\,\scriptstyle\mathfrak{Re}}$"
labels["BIm_x"] = r"${B_x}_{\,\scriptstyle\mathfrak{Im}}$"
labels["BIm_y"] = r"${B_y}_{\,\scriptstyle\mathfrak{Im}}$"
labels["BIm_z"] = r"${B_z}_{\,\scriptstyle\mathfrak{Im}}$"
labels["F_x"] = r"$\left<{F_x}\right>_{t}$"
labels["F_y"] = r"$\left<{F_y}\right>_{t}$"
labels["F_z"] = r"$\left<{F_z}\right>_{t}$"
# VReGrad_x VReGrad_y VReGrad_z VImGrad_x VImGrad_y VImGrad_z VRe VIm sigma mur

colors = dict()
colors["jRe_x"] = "hzdr-orange"
colors["jRe_y"] = "hzdr-yellow"
colors["jRe_z"] = "hzdr-green"
colors["jIm_x"] = "hzdr-blue"
colors["jIm_y"] = "hzdr-purple"
colors["jIm_z"] = "hzdr-red"
colors["BRe_x"] = "hzdr-orange"
colors["BRe_y"] = "hzdr-yellow"
colors["BRe_z"] = "hzdr-green"
colors["BIm_x"] = "hzdr-blue"
colors["BIm_y"] = "hzdr-purple"
colors["BIm_z"] = "hzdr-red"
colors["F_x"] = "hzdr-orange"
colors["F_y"] = "hzdr-yellow"
colors["F_z"] = "hzdr-green"
# VReGrad_x VReGrad_y VReGrad_z VImGrad_x VImGrad_y VImGrad_z VRe VIm sigma mur

markers = dict()

markers["jRe_x"] = "o"
markers["jRe_y"] = "v"
markers["jRe_z"] = "^"
markers["jIm_x"] = "s"
markers["jIm_y"] = "d"
markers["jIm_z"] = "*"
markers["BRe_x"] = "o"
markers["BRe_y"] = "v"
markers["BRe_z"] = "^"
markers["BIm_x"] = "s"
markers["BIm_y"] = "d"
markers["BIm_z"] = "*"
markers["F_x"] = "o"
markers["F_y"] = "v"
markers["F_z"] = "^"
# VReGrad_x VReGrad_y VReGrad_z VImGrad_x VImGrad_y VImGrad_z VRe VIm sigma mur

# --------------------------------------------------------------------------- #
# --- Test ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

def fig(case, opmesh, ofmesh, freq, line, name, fields,
        scaleX=1.0, scaleY=1.0):

    fig = plt.figure()

    def plot():

        axs = fig.add_subplot(111)

        axs.minorticks_on()

        axs.set_xlim([-75,75])
        #axs.set_ylim([0,30])

        axs.set_xlabel(alabels[line[:-1]])
        #axs.set_ylabel(labelAxisZ)

        #axs.set_aspect("equal")

        opData = data["Opera3D"][case][opmesh][freq][line]
        ofData = data["eddyCurrentFoam"][case][ofmesh][freq][line]

        opsx = scales["Opera3D"][line[:-1]]
        ofsx = scales["eddyCurrentFoam"][line[:-1]]

        for field in fields:

            opsy = scales["Opera3D"][field]
            ofsy = scales["eddyCurrentFoam"][field]

            color = colors[field]
            marker = markers[field]
            label = labels[field]

            axs.plot(opsx*scaleX*opData[line[:-1]], opsy*scaleY*opData[field],
                     color=color, linestyle="--")

            axs.plot(ofsx*scaleX*ofData[line[:-1]], ofsy*scaleY*ofData[field],
                     color=color, linestyle="-",
                     marker=marker, markevery=5, markersize=5,
                     markeredgecolor=color, markerfacecolor=color,
                     label=label)

        legendCols  = 3
        legendShift = 0.035 * (len(fields)/legendCols -1)

        axs.legend(bbox_to_anchor=(0.0, 1.05+legendShift, 1.0, 0.05+legendShift),
                  loc="upper center", ncol=legendCols,
                  mode="expand", borderaxespad=0.)

    plot()

    name = "plot_" + case + "_opm" + opmesh + "_ofm" + ofmesh + "_f" \
         + freq + "_line_" + line + "_" + name

    #fig.set_size_inches(sizeCompX, sizeCompY)
    fig.savefig(__dir__ + "/" + name + ".pdf", bbox_inches="tight")

    plt.close(fig)

fig("ortho", "1.000", "1.000", "1000", "y2", "j",
    ["jRe_x", "jRe_y", "jRe_z", "jIm_x", "jIm_y", "jIm_z"],
    scaleX=1e+3, scaleY=1e-6)

fig("ortho", "1.000", "1.000", "1000", "y2", "B",
    ["BRe_x", "BRe_y", "BRe_z", "BIm_x", "BIm_y", "BIm_z"],
    scaleX=1e+3, scaleY=1e+2)

fig("ortho", "1.000", "1.000", "1000", "y2", "F",
    ["F_x", "F_y", "F_z"],
    scaleX=1e+3, scaleY=1e-4)

fig("nonortho", "1.000", "0.125", "1000", "y2", "j",
    ["jRe_x", "jRe_y", "jRe_z", "jIm_x", "jIm_y", "jIm_z"],
    scaleX=1e+3, scaleY=1e-6)

fig("nonortho", "1.000", "0.125", "1000", "y2", "B",
    ["BRe_x", "BRe_y", "BRe_z", "BIm_x", "BIm_y", "BIm_z"],
    scaleX=1e+3, scaleY=1e+2)

fig("nonortho", "1.000", "0.125", "1000", "y2", "F",
    ["F_x", "F_y", "F_z"],
    scaleX=1e+3, scaleY=1e-4)

fig("nonortho", "1.000", "0.250", "1000", "y2", "j",
    ["jRe_x", "jRe_y", "jRe_z", "jIm_x", "jIm_y", "jIm_z"],
    scaleX=1e+3, scaleY=1e-6)

fig("nonortho", "1.000", "0.250", "1000", "y2", "B",
    ["BRe_x", "BRe_y", "BRe_z", "BIm_x", "BIm_y", "BIm_z"],
    scaleX=1e+3, scaleY=1e+2)

fig("nonortho", "1.000", "0.250", "1000", "y2", "F",
    ["F_x", "F_y", "F_z"],
    scaleX=1e+3, scaleY=1e-4)

fig("nonortho", "1.000", "0.375", "1000", "y2", "j",
    ["jRe_x", "jRe_y", "jRe_z", "jIm_x", "jIm_y", "jIm_z"],
    scaleX=1e+3, scaleY=1e-6)

fig("nonortho", "1.000", "0.375", "1000", "y2", "B",
    ["BRe_x", "BRe_y", "BRe_z", "BIm_x", "BIm_y", "BIm_z"],
    scaleX=1e+3, scaleY=1e+2)

fig("nonortho", "1.000", "0.375", "1000", "y2", "F",
    ["F_x", "F_y", "F_z"],
    scaleX=1e+3, scaleY=1e-4)

fig("nonortho", "1.000", "0.500", "1000", "y2", "j",
    ["jRe_x", "jRe_y", "jRe_z", "jIm_x", "jIm_y", "jIm_z"],
    scaleX=1e+3, scaleY=1e-6)

fig("nonortho", "1.000", "0.500", "1000", "y2", "B",
    ["BRe_x", "BRe_y", "BRe_z", "BIm_x", "BIm_y", "BIm_z"],
    scaleX=1e+3, scaleY=1e+2)

fig("nonortho", "1.000", "0.500", "1000", "y2", "F",
    ["F_x", "F_y", "F_z"],
    scaleX=1e+3, scaleY=1e-4)

fig("nonortho", "1.000", "0.750", "1000", "y2", "j",
    ["jRe_x", "jRe_y", "jRe_z", "jIm_x", "jIm_y", "jIm_z"],
    scaleX=1e+3, scaleY=1e-6)

fig("nonortho", "1.000", "0.750", "1000", "y2", "B",
    ["BRe_x", "BRe_y", "BRe_z", "BIm_x", "BIm_y", "BIm_z"],
    scaleX=1e+3, scaleY=1e+2)

fig("nonortho", "1.000", "0.750", "1000", "y2", "F",
    ["F_x", "F_y", "F_z"],
    scaleX=1e+3, scaleY=1e-4)

fig("nonortho", "1.000", "1.000", "1000", "y2", "j",
    ["jRe_x", "jRe_y", "jRe_z", "jIm_x", "jIm_y", "jIm_z"],
    scaleX=1e+3, scaleY=1e-6)

fig("nonortho", "1.000", "1.000", "1000", "y2", "B",
    ["BRe_x", "BRe_y", "BRe_z", "BIm_x", "BIm_y", "BIm_z"],
    scaleX=1e+3, scaleY=1e+2)

fig("nonortho", "1.000", "1.000", "1000", "y2", "F",
    ["F_x", "F_y", "F_z"],
    scaleX=1e+3, scaleY=1e-4)

fig("nonortho", "1.000", "1.500", "1000", "y2", "j",
    ["jRe_x", "jRe_y", "jRe_z", "jIm_x", "jIm_y", "jIm_z"],
    scaleX=1e+3, scaleY=1e-6)

fig("nonortho", "1.000", "1.500", "1000", "y2", "B",
    ["BRe_x", "BRe_y", "BRe_z", "BIm_x", "BIm_y", "BIm_z"],
    scaleX=1e+3, scaleY=1e+2)

fig("nonortho", "1.000", "1.500", "1000", "y2", "F",
    ["F_x", "F_y", "F_z"],
    scaleX=1e+3, scaleY=1e-4)

fig("nonortho", "1.000", "2.000", "1000", "y2", "j",
    ["jRe_x", "jRe_y", "jRe_z", "jIm_x", "jIm_y", "jIm_z"],
    scaleX=1e+3, scaleY=1e-6)

fig("nonortho", "1.000", "2.000", "1000", "y2", "B",
    ["BRe_x", "BRe_y", "BRe_z", "BIm_x", "BIm_y", "BIm_z"],
    scaleX=1e+3, scaleY=1e+2)

fig("nonortho", "1.000", "2.000", "1000", "y2", "F",
    ["F_x", "F_y", "F_z"],
    scaleX=1e+3, scaleY=1e-4)

fig("nonortho", "1.000", "2.500", "1000", "y2", "j",
    ["jRe_x", "jRe_y", "jRe_z", "jIm_x", "jIm_y", "jIm_z"],
    scaleX=1e+3, scaleY=1e-6)

fig("nonortho", "1.000", "2.500", "1000", "y2", "B",
    ["BRe_x", "BRe_y", "BRe_z", "BIm_x", "BIm_y", "BIm_z"],
    scaleX=1e+3, scaleY=1e+2)

fig("nonortho", "1.000", "2.500", "1000", "y2", "F",
    ["F_x", "F_y", "F_z"],
    scaleX=1e+3, scaleY=1e-4)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #