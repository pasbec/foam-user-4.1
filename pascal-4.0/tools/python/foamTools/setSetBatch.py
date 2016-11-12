#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Module template
# March 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

from .ioInfo import ioBase

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

# --------------------------------------------------------------------------- #
# --- Class definitions ----------------------------------------------------- #
# --------------------------------------------------------------------------- #

class setSetBatch(object):

    def __init__(self, fileName=None):

        self.nameStringMaxLenth = 60

        self.io = ioBase(fileName)

    # ----------------------------------------------------------------------- #

    def _topoSetString(self, string):

        return string.ljust(10, " ")

    # ----------------------------------------------------------------------- #

    def _nameString(self, string):

        return string.ljust(self.nameStringMaxLenth, " ")

    # ----------------------------------------------------------------------- #

    def _actionString(self, string):

        return string.ljust(7, " ")

    # ----------------------------------------------------------------------- #

    def _writeTopoSetString(self, string):

        self.io.write(self._topoSetString(string), end=" ")

    # ----------------------------------------------------------------------- #

    def _writeNameString(self, string):

        self.io.write(self._nameString(string), end=" ")

    # ----------------------------------------------------------------------- #

    def _writeActionString(self, string):

        self.io.write(self._actionString(string), end=" ")

    # ----------------------------------------------------------------------- #

    def _doManual(self, topo, name, action, source):

# TODO: Check for valid actions?
# TODO: Check for valid sources?

        self._writeTopoSetString(topo + "Set")
        self._writeNameString(topo + "Set" + "_" + name)
        self._writeActionString(action)
        self.io.write(source)

    # ----------------------------------------------------------------------- #

    def _doFromStl(self, topo, name, stl, tol=2e-4, i=False, o=False):

        # Define write for boolean-strings for OpenFOAM
        def bstr(b): return "true" if b else "false"

# TODO: Check name for NO spaces

        self._writeTopoSetString("pointSet")
        self._writeNameString("pointSet" + "_" + name)
        self._writeActionString("new")
        self.io.write("surfaceToPoint" + " ", end=" ")
        self.io.write("\"" + stl + "\"" + str(tol) + " " + bstr(i) + " " + bstr(o))

        self._writeTopoSetString(topo + "Set")
        self._writeNameString(topo + "Set" + "_" + name)
        self._writeActionString("new")
        self.io.write("pointTo" + topo.capitalize() + " ", end=" ")
        self.io.write("pointSet" + "_" + name + "all")

        self._writeTopoSetString("pointSet")
        self._writeNameString("pointSet" + "_" + name)
        self._writeActionString("remove")
        self.io.write()

    # ----------------------------------------------------------------------- #

    _do = {"manual": _doManual,
           "fromStl": _doFromStl}

    # ----------------------------------------------------------------------- #

    def cellSet(self, do, name, *args, **kwargs):

# TODO: Check name for NO spaces

        self._do[do](self, "cell", name, *args, **kwargs)

    # ----------------------------------------------------------------------- #

    def manual(self, string, end="\n"):

        self.io.write(string, end=end)

# --------------------------------------------------------------------------- #
# --- End of module --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

