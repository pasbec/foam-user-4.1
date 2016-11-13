#!/usr/bin/python
#
# October 2016
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

__name__
__path__ = os.path.realpath(__file__)
__base__ = os.path.basename(__path__)
__dir__ = os.path.dirname(__path__)
__head__ = os.path.splitext(__base__)[0]

import matplotlib as mpl
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------- #
# --- Function -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def latexify(
        fig_width=None, fig_height=None, columns=1,
        fontsize=8, fontfamily="serif"
    ):
    """Set up matplotlib"s RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, mm
    fig_height : float,  optional, mm
    columns : int {1, 2}, optional
    fontsize : int, optional
    fontfamily : str, optional
    """

    #inch = 25.4

    #if fig_width is not None:
        #fig_width = fig_width / inch

    #if fig_height is not None:
        #fig_height = fig_width / inch

    ## Code adapted from http://www.scipy.org/Cookbook/Matplotlib/LaTeX_Examples

    ## Width and max height in inches for IEEE journals taken from
    ## computer.org/cms/Computer.org/Journal%20templates/transactions_art_guide.pdf

    #assert(columns in [1,2])

    #if fig_width is None:
        #fig_width = 3.39 if columns==1 else 6.9 # width in inches

    #if fig_height is None:
        #golden_mean = ((5**0.5)-1.0)/2.0    # Aesthetic ratio
        #fig_height = fig_width*golden_mean # height in inches

    #MAX_HEIGHT_INCHES = 8.0
    #if fig_height > MAX_HEIGHT_INCHES:
        #print("WARNING: fig_height too large:" + fig_height +
              #"so will reduce to" + MAX_HEIGHT_INCHES + "inches.")
        #fig_height = MAX_HEIGHT_INCHES

    ## Adjust size

    #mpl.rcParams["figure.figsize"] = [fig_width,fig_height]

    # Adjust settings

    mpl.rcParams["backend"] = "ps"

    mpl.rcParams["text.usetex"] = True
    mpl.rcParams["text.latex.unicode"] = True
    #mpl.rcParams["text.latex.preamble"] = ["\usepackage{gensymb}"]

    mpl.rcParams["font.size"] = fontsize
    mpl.rcParams["font.family"] = fontfamily

    mpl.rcParams["axes.labelsize"] = fontsize
    mpl.rcParams["axes.titlesize"] = fontsize

    mpl.rcParams["legend.fontsize"] = fontsize

    mpl.rcParams["xtick.labelsize"] = fontsize
    mpl.rcParams["ytick.labelsize"] = fontsize

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

