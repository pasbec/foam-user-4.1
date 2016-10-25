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

# pyFoam modules
from PyFoam.RunDictionary.ParsedParameterFile import ParsedParameterFile
from PyFoam.Basics.DataStructures import Vector,Field,Dimension,DictProxy,TupleProxy,Tensor,SymmTensor,Unparsed,UnparsedList,Codestream,DictRedirection,BinaryBlob,BinaryList,BoolProxy

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

def tryKey(d, k):

    # return : boolean for key status
    #
    # d      : dictionary
    # k      : key of dictionary to try

    # TODO [Low]: Improve descriptions

    try:
        assert d[k]==d[k]

    except:
        return False

    return True



def tryMandatory(d, k):

    # return : boolean for key status
    #
    # d      : dictionary
    # k      : mandatory key of dictionary

    # TODO [Low]: Improve descriptions

    if tryKey(d, k):
        return d[k]

    else:
        raise KeyError('mandatory key "' + k + '" is missing')



def tryOrDefault(d, k, vDef):

    # return : boolean for key status
    #
    # d      : dictionary
    # k      : key of dictionary to try
    # vDef   : default value of key

    # TODO [Low]: Improve descriptions

    if tryKey(d, k):
        return d[k]

    else:
        return vDef



def readStr(v):

    # return : stripped dictionary value (without '"')
    #
    # v      : dictionary key value of type str

    # TODO [Low]: Improve descriptions

    if type(v) == str:
        return v.strip('"')
    else:
        try:
            return str(v)
        except:
            raise TypeError('key value ' + str(v) + ' is no string')



def readInt(v):

    # return : integer dictionary value
    #
    # v      : dictionary key value of type int

    # TODO [Low]: Improve descriptions

    if type(v) == int:
        return v

    if type(v) == str:
        return int(readStr(v))

    else:
      raise TypeError('key value ' + str(v) + ' is not an integer')



def readFloat(v):

    # return : float dictionary value
    #
    # v      : dictionary key value of type float

    # TODO [Low]: Improve descriptions

    if type(v) == float:
        return v

    if type(v) == int:
        return float(v)

    if type(v) == str:
        return float(readStr(v))

    else:
      raise TypeError('key value ' + str(v) + ' is not an float')



def readDim(v):

    # return : integer dictionary value
    #
    # v      : dictionary key value of type Dimension

    # TODO [Low]: Improve descriptions

    if type(v) == Dimension:
        return v

    if type(v) == str:
        dim = readStr(v).strip('[]')
        dims = dim.split()
        print dims
        return Dimension(*dims)

    else:
      raise TypeError('key value ' + str(v) + ' has no dimensional type')



def readBool(v):

    # return : interpreted dictionary boolean value
    #
    # v      : dictionary key value of type BoolProxy

    # TODO [Low]: Improve descriptions

    if type(v) == BoolProxy:
        return v

    elif type(v) == str:
        trueKeys = ['On', 'on', 'True', 'true', 'Yes', 'yes']
        falseKeys  = ['Off', 'off', 'False', 'false', 'No', 'no']
        for k in trueKeys:
            if readStr(v) == k: return BoolProxy(val=True)
        for k in falseKeys:
            if readStr(v) == k: return BoolProxy(val=False)
        raise TypeError('keyValue ' + str(v) + ' is no boolean')

    elif type(v) == int:
        if v == 1: return BoolProxy(val=True)
        elif v == 0: return BoolProxy(val=False)

    else:
      raise TypeError('key value ' + str(v) + ' is no boolean')



def readSubDict(v):

    # return : integer dictionary value
    #
    # v      : dictionary key value of type DictProxy

    # TODO [Low]: Improve descriptions

    if type(v) == DictProxy:
        return v

    else:
      raise TypeError('key value ' + str(v) + ' is no (sub-)dictionary')

# --------------------------------------------------------------------------- #
# --- End of module --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

