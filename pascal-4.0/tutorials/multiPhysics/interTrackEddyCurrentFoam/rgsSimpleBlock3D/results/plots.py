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

def readData(set, cases, frequencies, lines, meshes):

    setData = dict()
    setFields = dict()
    setMeshErr = dict()
    setErrNorms = dict()
    setErrNormMeshes = dict()

    for case in cases:

        setData[case] = dict()
        setFields[case] = dict()
        setMeshErr[case] = dict()
        setErrNorms[case] = dict()
        setErrNormMeshes[case] = dict()

        for freq in frequencies:

            setData[case][freq] = dict()
            setFields[case][freq] = dict()
            setMeshErr[case][freq] = dict()
            setErrNorms[case][freq] = dict()
            setErrNormMeshes[case][freq] = dict()

            for line in lines:

                setData[case][freq][line] = dict()
                setFields[case][freq][line] = dict()
                setMeshErr[case][freq][line] = dict()
                setErrNorms[case][freq][line] = dict()
                setErrNormMeshes[case][freq][line] = dict()

                meshRef = False

                for mesh in reversed(sorted(meshes)):

                    fileName = __dir__ + "/data_" + set \
                             + "_" + case + "_f" + freq + "_line_" + line \
                             + "_m" + mesh + ".dat"

                    if os.path.isfile(fileName):

                        if not meshRef: meshRef = mesh

                        print("Reading file: " + fileName)

                        setData[case][freq][line][mesh] = \
                            np.genfromtxt(fileName, comments='#', names=True)

                        data = setData[case][freq][line][mesh]
                        dataRef = setData[case][freq][line][meshRef]

                        setFields[case][freq][line][mesh] = \
                            data.dtype.names

                        setMeshErr[case][freq][line][mesh] = \
                            np.zeros(data.shape, data.dtype)


                        for field in setFields[case][freq][line][mesh]:

                            setMeshErr[case][freq][line][mesh][field] = \
                                np.abs(data[field] - dataRef[field])

                            if not field in setErrNorms[case][freq][line]:

                                setErrNorms[case][freq][line][field] = dict()
                                setErrNormMeshes[case][freq][line][field] = []

                            error = setMeshErr[case][freq][line][mesh][field]

                            errNorms = setErrNorms[case][freq][line][field]

                            errNormMeshes = setErrNormMeshes[case][freq][line][field]
                            meshStored = False

                            for norm in norms:

                                if not norm in errNorms:

                                    errNorms[norm] = []

                                else:

                                    errNorms[norm].insert(0, norms[norm](dataRef[field], error))

                                    if not meshStored:

                                        errNormMeshes.insert(0, mesh)
                                        meshStored=True

    return setData, setFields, setMeshErr, setErrNorms, setErrNormMeshes

# --------------------------------------------------------------------------- #
# --- Data ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

cases = ["ortho", "nonortho"]
frequencies = ["1000", "10000", "100000"]
lines = ["x1", "y1", "y2", "z1"]
meshes = ["0.125", "0.250", "0.375", "0.500", "0.750", "1.000",
          "1.500", "2.000", "2.500"]

# --------------------------------------------------------------------------- #

data = dict()
fields = dict()
error = dict()
errNorms = dict()
errNormMeshes = dict()
scales = dict()

# --------------------------------------------------------------------------- #

set = "Opera3D"

data[set], fields[set], error[set], errNorms[set], errNormMeshes[set] = \
    readData(set, cases, frequencies, lines, meshes)

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

data[set], fields[set], error[set], errNorms[set], errNormMeshes[set] = \
    readData(set, cases, frequencies, lines, meshes)

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
scales[set]["VReGrad_x"] = 1.0
scales[set]["VReGrad_y"] = 1.0
scales[set]["VReGrad_z"] = 1.0
scales[set]["VImGrad_x"] = 1.0
scales[set]["VImGrad_y"] = 1.0
scales[set]["VImGrad_z"] = 1.0
scales[set]["VRe"] = 1.0
scales[set]["VIm"] = 1.0
scales[set]["sigma"] = 1.0
scales[set]["mur"] = 1.0

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
labels["VReGrad_x"] = r"${\nabla \phi_x}_{\,\scriptstyle\mathfrak{Re}}$"
labels["VReGrad_y"] = r"${\nabla \phi_y}_{\,\scriptstyle\mathfrak{Re}}$"
labels["VReGrad_z"] = r"${\nabla \phi_z}_{\,\scriptstyle\mathfrak{Re}}$"
labels["VImGrad_x"] = r"${\nabla \phi_x}_{\,\scriptstyle\mathfrak{Im}}$"
labels["VImGrad_y"] = r"${\nabla \phi_y}_{\,\scriptstyle\mathfrak{Im}}$"
labels["VImGrad_z"] = r"${\nabla \phi_z}_{\,\scriptstyle\mathfrak{Im}}$"
labels["VRe"] = r"${\phi}_{\,\scriptstyle\mathfrak{Re}}$"
labels["VIm"] = r"${\phi}_{\,\scriptstyle\mathfrak{Im}}$"
labels["sigma"] = r"${\sigma}$"
labels["mur"] = r"${\mu_\mathrm{r}}$"

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
colors["F_y"] = "hzdr-blue"
colors["F_z"] = "hzdr-green"
colors["VReGrad_x"] = "hzdr-orange"
colors["VReGrad_y"] = "hzdr-yellow"
colors["VReGrad_z"] = "hzdr-green"
colors["VImGrad_x"] = "hzdr-blue"
colors["VImGrad_y"] = "hzdr-purple"
colors["VImGrad_z"] = "hzdr-red"
colors["VRe"] = "hzdr-orange"
colors["VIm"] = "hzdr-blue"
colors["sigma"] = "hzdr-blue"
colors["mur"] = "hzdr-blue"

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
markers["VReGrad_x"] = "o"
markers["VReGrad_y"] = "v"
markers["VReGrad_z"] = "^"
markers["VImGrad_x"] = "s"
markers["VImGrad_y"] = "d"
markers["VImGrad_z"] = "*"
markers["VRe"] = "o"
markers["VIm"] = "s"
markers["sigma"] = "o"
markers["mur"] = "o"

# --------------------------------------------------------------------------- #
# --- Test ------------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

def fig(case, freq, line, mesh, name, fields, op=True,
        scaleX=1.0, scaleY=1.0, shiftLegend=0.0):

    fig = plt.figure()

    opmesh = "1.000"

    def plot():

        axs = fig.add_subplot(111)

        axs.minorticks_on()

        axs.set_xlim([-75,75])
        #axs.set_ylim([0,30])

        axs.set_xlabel(alabels[line[:-1]])
        #axs.set_ylabel(labelAxisZ)

        #axs.set_aspect("equal")

        if op:

            opsx = scales["Opera3D"][line[:-1]]
            opData = data["Opera3D"][case][freq][line][opmesh]

        ofData = data["eddyCurrentFoam"][case][freq][line][mesh]
        ofsx = scales["eddyCurrentFoam"][line[:-1]]

        for field in fields:

            color = colors[field]
            marker = markers[field]
            label = labels[field]

            if op:

                opsy = scales["Opera3D"][field]

                axs.plot(opsx*scaleX*opData[line[:-1]],
                         opsy*scaleY*opData[field],
                         color=color, linestyle="--")

            ofsy = scales["eddyCurrentFoam"][field]

            axs.plot(ofsx*scaleX*ofData[line[:-1]],
                     ofsy*scaleY*ofData[field],
                     color=color, linestyle="-",
                     marker=marker, markevery=5, markersize=5,
                     markeredgecolor=color, markerfacecolor=color,
                     label=label)

        legendCols  = 3
        legendPos = shiftLegend + 0.035 * (len(fields)/legendCols - 1)

        axs.legend(bbox_to_anchor=(0.0, 1.05+legendPos, 1.0, 0.05+legendPos),
                   loc="upper center", ncol=legendCols,
                   mode="expand", borderaxespad=0.)

    plot()

    fileName = "plot_" + case + "_f" + freq + "_line_" + line + "_m" + mesh

    if op: fileName += "-op" + opmesh

    fileName += "_" + name

    #fig.set_size_inches(sizeCompX, sizeCompY)
    fig.savefig(__dir__ + "/" + fileName + ".pdf", bbox_inches="tight")

    plt.close(fig)

#for mesh in ["1.000"]:

    #fig("ortho", "1000", "y2", mesh, "j",
        #["jRe_x", "jRe_y", "jRe_z", "jIm_x", "jIm_y", "jIm_z"],
        #scaleX=1e+3, scaleY=1e-6)

    #fig("ortho", "1000", "y2", mesh, "B",
        #["BRe_x", "BRe_y", "BRe_z", "BIm_x", "BIm_y", "BIm_z"],
        #scaleX=1e+3, scaleY=1e+2)

    #fig("ortho", "1000", "y2", mesh, "F",
        #["F_x", "F_y", "F_z"],
        #scaleX=1e+3, scaleY=1e-4)

    #fig("ortho", "1000", "y2", mesh, "V",
        #["VRe", "VIm"], op=False,
        #scaleX=1e+3, scaleY=1.0, shiftLegend=0.01)

    #fig("ortho", "1000", "y2", mesh, "VGrad",
        #["VReGrad_x", "VReGrad_y", "VReGrad_z",
         #"VImGrad_x", "VImGrad_y", "VImGrad_z"], op=False,
        #scaleX=1e+3, scaleY=11.0)

    #fig("ortho", "1000", "y2", mesh, "sigma",
        #["sigma"], op=False,
        #scaleX=1e+3, scaleY=1.0, shiftLegend=0.02)

    #fig("ortho", "1000", "y2", mesh, "mur",
        #["mur"], op=False,
        #scaleX=1e+3, scaleY=1.0, shiftLegend=0.01)

#for mesh in ["1.000"]:

    #fig("nonortho", "1000", "y2", mesh, "j",
        #["jRe_x", "jRe_y", "jRe_z", "jIm_x", "jIm_y", "jIm_z"],
        #scaleX=1e+3, scaleY=1e-6)

    #fig("nonortho", "1000", "y2", mesh, "B",
        #["BRe_x", "BRe_y", "BRe_z", "BIm_x", "BIm_y", "BIm_z"],
        #scaleX=1e+3, scaleY=1e+2)

    #fig("nonortho", "1000", "y2", mesh, "F",
        #["F_x", "F_y", "F_z"],
        #scaleX=1e+3, scaleY=1e-4)

    #fig("nonortho", "1000", "y2", mesh, "V",
        #["VRe", "VIm"], op=False,
        #scaleX=1e+3, scaleY=1.0, shiftLegend=0.01)

    #fig("nonortho", "1000", "y2", mesh, "VGrad",
        #["VReGrad_x", "VReGrad_y", "VReGrad_z",
         #"VImGrad_x", "VImGrad_y", "VImGrad_z"], op=False,
        #scaleX=1e+3, scaleY=11.0)

    #fig("nonortho", "1000", "y2", mesh, "sigma",
        #["sigma"], op=False,
        #scaleX=1e+3, scaleY=1.0, shiftLegend=0.02)

    #fig("nonortho", "1000", "y2", mesh, "mur",
        #["mur"], op=False,
        #scaleX=1e+3, scaleY=1.0, shiftLegend=0.02)

# --------------------------------------------------------------------------- #

def test(case, freq, line, mesh, field,
        scaleX=1.0, scaleY=1.0, shiftLegend=0.0):

    fig = plt.figure()

    def plot():

        axs = fig.add_subplot(111)

        axs.minorticks_on()

        axs.set_xscale("log")
        axs.set_yscale("log")

        data = errNorms["eddyCurrentFoam"][case][freq][line][field]["1"]
        print(data)
        meshes = errNormMeshes["eddyCurrentFoam"][case][freq][line][field]
        print(meshes)

        axs.plot(data, meshes)

    plot()

    fileName = "test"

    #fig.set_size_inches(sizeCompX, sizeCompY)
    fig.savefig(__dir__ + "/" + fileName + ".pdf", bbox_inches="tight")

    plt.close(fig)

test("nonortho", "1000", "y2", "1.000", "F_y",
    scaleX=1e+3, scaleY=1.0, shiftLegend=0.02)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
