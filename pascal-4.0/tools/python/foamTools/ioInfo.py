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

def fileGetScriptPath():

    return os.path.dirname(os.path.realpath(sys.argv[0]))



def fileGetPath(fileName):

    if os.path.isabs(fileName):

        return fileName

    else:

        return fileGetScriptPath() + "/" + fileName



def objectIndent(cString, iLevel=0, iChar=' ', iCount=4):

    # return  : indented content string
    #
    # cString : content string to indent
    # iLevel  : indentation level
    # iChar   : indentation character
    # iCount  : indentation character count

    return iLevel * iChar * iCount + cString



def objectHeader(name, cl, time=""):

    # return : header string
    #
    # time   : time instant
    # name   : name of foam object
    # cl     : class of foam object (e.g. volScalarField)

    # Define short indented line with line break
    def i(iL,cS,eS='\n'): return objectIndent(cS + eS,iLevel=iL)

    # Assemble header string
    r  = ''
    r += i(0, '/*--------------------------------*- C++ -*----------------------------------*\\')
    r += i(0, '| =========                 |                                                 |')
    r += i(0, '| \\\\      /  F ield         | foam-extend: Open Source CFD                    |')
    r += i(0, '|  \\\\    /   O peration     | Version:     4.0                                |')
    r += i(0, '|   \\\\  /    A nd           | Web:         http://www.extend-project.de       |')
    r += i(0, '|    \\\\/     M anipulation  | For copyright notice see file Copyright         |')
    r += i(0, '\*---------------------------------------------------------------------------*/')
    r += i(0, 'FoamFile')
    r += i(0, '{')
    r += i(1, 'version     2.0;')
    r += i(1, 'format      ascii;')
    r += i(1, 'class       ' + cl + ';')
    if time: r += i(1, 'location    "' + time + '"' + ';')
    r += i(1, 'object      ' + name + ';')
    r += i(0, '}')
    r += i(0, '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //')
    r += i(0, '')

    return r



def objectFooter():

    # return : footer string

    # Define short indented line with line break
    def i(iL,cS,eS='\n'): return objectIndent(cS + eS,iLevel=iL)

    # Assemble footer string
    r  = ''
    r += i(0, '// ************************************************************************* //')
    r += i(0, '')

    return r

# --------------------------------------------------------------------------- #
# --- End of module --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

