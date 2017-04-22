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

dataBaseName = "data"
plotBaseName = "plot"

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

set = "eddyCurrentFoam"

data[set] = dict()
names[set] = dict()

meshes = ["1.000"]
lines = ["x1", "y1", "y2", "z1"]

# Read data
for mesh in meshes:

    data[set][mesh] = dict()
    names[set][mesh] = dict()

    for line in lines:

        fileName = __dir__ + "/" + dataBaseName + "_" + set \
                 + "_" + mesh + "_line_" + line + ".dat"

        data[set][mesh][line] = np.genfromtxt(fileName, comments='#',
                                              names=True, excludelist=["s2"])
        names[set][mesh][line] = data[set][mesh][line].dtype.names

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

        ax.plot(data["eddyCurrentFoam"]["1.000"]["y2"]["s"],
                data["eddyCurrentFoam"]["1.000"]["y2"]["VRe"],
                label="VRe")

        #ax.set_xlim([0,30])
        #ax.set_ylim([0,30])

        #ax.set_xlabel(labelAxisR)
        #ax.set_ylabel(labelAxisZ)

        #ax.set_aspect("equal")

        #set = "Analytical"
        #c = ax.contour(R[set], Z[set], F[set][1],
                       #levels=levels[1])
        #cl = ax.clabel(c, c.levels[0::2],
                       #inline=True, fmt="%g", fontsize=fontsize)

        #[l.set_bbox(dict(facecolor="white", edgecolor="none", pad=2)) for l in cl]
        #[l.set_text('{:n}'.format(float(l.get_text()))) for l in cl]

    ax(fig, axs, "test")

    #fig.set_size_inches(sizeCompX, sizeCompY)
    #fig.savefig(__dir__+"/"+plotBaseName+name+".pdf", bbox_inches="tight")
    fig.savefig(__dir__+ "/" + name + ".pdf", bbox_inches="tight")

fig(plots, "test")

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
