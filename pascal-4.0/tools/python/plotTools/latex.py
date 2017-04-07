#!/usr/bin/python
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

import locale as lc
import matplotlib as mpl
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------- #
# --- Function -------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

def latexify(
        fig_width=None, fig_height=None, dpi=None,
        fontsize=None, fontfamily=None, locale=None
    ):
    """Set up matplotlib"s RC params for LaTeX plotting.
    Call this before plotting a figure.

    Parameters
    ----------
    fig_width : float, optional, inch
    fig_height : float,  optional, inch
    dpi : int, optional,
    fontsize : int, optional
    fontfamily : str, optional
    """

    inch = 25.4

    ## Adjust size and dpi
    if fig_width is not None or fig_height is not None:

        mpl.rcParams["figure.figsize"] = [fig_width,fig_height]

    if dpi is not None:

        mpl.rcParams["figure.dpi"] = dpi
        mpl.rcParams["savefig.dpi"] = dpi

    # Adjust settings

    mpl.rcParams["backend"] = "ps"

    mpl.rcParams["text.usetex"] = True
    mpl.rcParams["text.latex.unicode"] = True
    #mpl.rcParams["text.latex.preamble"] = ["\usepackage{gensymb}"]

    if fontfamily is not None:

        mpl.rcParams["font.family"] = fontfamily

    if fontsize is not None:

        mpl.rcParams["font.size"] = fontsize

        mpl.rcParams["axes.labelsize"] = fontsize
        mpl.rcParams["axes.titlesize"] = fontsize

        mpl.rcParams["legend.fontsize"] = fontsize

        mpl.rcParams["xtick.labelsize"] = fontsize
        mpl.rcParams["ytick.labelsize"] = fontsize

    # Adjust locale
    if locale is not None:

        lc.setlocale(lc.LC_ALL, locale)
        mpl.rcParams['axes.formatter.use_locale'] = True

    # Adjust savefig?


# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

