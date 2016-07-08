#!/usr/bin/python
# 
# March 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# TODO [High]: Allow other primitive patch types than 'fixedValue'

# TODO [Low]: Rework function descriptions

# TODO [Low]: Reorganize some functions in modules

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

import re

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

def fieldClassToAtomicTypeSize(cl):

    # return : atomic data type and its size of foam object (e.g. scalar)
    #
    # cl     : class of foam object (e.g. volScalarField)

    # Predefined mappings
    scalarTypes = ['scalarField', 'volScalarField', 'surfaceScalarField', 'pointScalarField', 'tetPointScalarField']
    vectorTypes = ['vectorField', 'volVectorField', 'surfaceVectorField', 'pointVectorField', 'tetPointVectorField']

    if any(cl in t for t in scalarTypes): return 'scalar', 1
    elif any(cl in t for t in vectorTypes): return 'vector', 3



def patchBaseTypeToPrimitiveValued(bT):

    # return : primitive type and boolean indicating necessary
    #          patch fixedValue data
    #
    # bT     : base type of patch

    # Predefined mappings
    valueTypes = [('patch', 'fixedValue'), ('wall', 'fixedValue'), ('processor', 'fixedValue')]
    passTypes  = ['empty', 'cyclic', 'symmetryPlane', 'wedge']

    for t in valueTypes:
        if bT == t[0]: return t[1], True
    for t in passTypes:
        if bT == t: return t, False

# --------------------------------------------------------------------------- #
# --- Argument parsing ------------------------------------------------------ #
# --------------------------------------------------------------------------- #



# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #



# --------------------------------------------------------------------------- #
# --- Main program sequence ------------------------------------------------- #
# --------------------------------------------------------------------------- #



# --------------------------------------------------------------------------- #
# --- End of module --------------------------------------------------------- #
# --------------------------------------------------------------------------- #
