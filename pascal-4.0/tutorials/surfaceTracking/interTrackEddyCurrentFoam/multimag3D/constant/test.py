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

sys.path.append(os.environ['FOAM_USER_TOOLS']+'/python')

import math as m
import numpy as np

import matplotlib.pyplot as plt

import plotTools.latex as latex
import plotTools.hzdr as hzdr

from foamTools.expansion import expansion_n_ds_root, expansion_n_ds

import scipy as sp
import scipy.optimize as spo

import time

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# --- Plot settings --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

latex.latexify()
hzdr.colors()

# --------------------------------------------------------------------------- #

e = 2.0
n = 100
l = 120.0

ds = 1.0

n = expansion_n_ds(e, ds, l=l)

print(n)

x = np.linspace(2.0,100.0,100)
y = expansion_n_ds_root(e, ds, x, l=l)

plt.plot(x,y)

plt.show(block=False)
time.sleep(30)

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
