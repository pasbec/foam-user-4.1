#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Module with helper functions to write a blockdict with python
# July 2015
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import os, sys

import math as m
import numpy as np
import scipy as sp
import scipy.optimize as spo

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

def expansion_n_ds_root(e, ds, n, l=None):

    ls = ds

    if l: ls /= l

    f = pow(e, 1.0/(n - 1))

    return ls*(1.0 - pow(f, n)) - (1.0 - pow(f, 1))



def expansion_n_ds(e, ds, l=None, tol=1e-16, maxiter=100):

    ls = ds

    if l: ls /= l

    def N(n): return expansion_n_ds_root(e=e, ds=ls, n=n)

    n0 = 2

    return int(spo.newton(func=N, x0=n0, tol=tol, maxiter=maxiter))



def expansion_f_ds_root(n, ds, f, l=None):

    ls = ds

    if l: ls /= l

    return ls*(1.0 - pow(f, n)) - (1.0 - pow(f, 1))



def expansion_f_ds(n, ds, l=None, tol=1e-16, maxiter=100):

    if (type(n) == int) and (n > 1):

        ls = ds

        if l: ls /= l

        def F(f): return expansion_f_ds_root(n=n, ds=ls, f=f)

        # Make sure we start root search on the right side of the local max
        f0 = 1.0 + 2.0*(pow(1.0/(ls*n), 1.0/(n - 1)) - 1.0)

        return spo.newton(func=F, x0=f0, tol=tol, maxiter=maxiter)

    elif (type(n) == int) and (n == 1):

        return 1.0

    else:

        raise ValueError("n must be an integer, which is larger than 0")



def expansion_f_e(n, e):

    if (type(n) == int) and (n > 1):

        return pow(e, 1.0/(n - 1))

    elif (type(n) == int) and (n == 1):

        return 1.0

    else:

        raise ValueError("n must be an integer, which is larger than 0")



def expansion_e_ds(n, ds, l=None):

    if (type(n) == int) and (n > 1):

        ls = ds

        if l: ls /= l

        f = expansion_f_ds(n, ls)

        return pow(f, n - 1)

    elif (type(n) == int) and (n == 1):

        return 1.0

    else:

        raise ValueError("n must be an integer, which is larger than 0")



def expansion_e_f(n, f):

    if (type(n) == int) and (n > 1):

        return pow(f, n - 1)

    elif (type(n) == int) and (n == 1):

        return 1.0

    else:

        raise ValueError("n must be an integer, which is larger than 0")



def expansion_ds_e(n, e, l=None):

    if (type(n) == int) and (n > 1):

        f = expansion_f_e(n, e)

        if f == 1.0:

            ds = 1.0/n

        else:

            ds = (1.0 - pow(f, 1))/(1.0 - pow(f, n))

        if l: ds *= l

        return ds

    else:

        raise ValueError("n must be an integer, which is larger than 1")



def expansion_de_e(n, e, l=None):

    if (type(n) == int) and (n > 1):

        f = expansion_f_e(n, e)

        if f == 1.0:

            de = 1.0/n

        else:

            de = 1.0 - (1.0 - pow(f, n - 1))/(1.0 - pow(f, n))

        if l: de *= l

        return de

    else:

        raise ValueError("n must be an integer, which is larger than 1")



def expansion(n, e, l=None):

    if (type(n) == int) and (n > 0):

        f = expansion_f_e(n, e)

        t = np.zeros(n+1)

        if f == 1.0:

            for i in range(n+1): t[i] = float(i)/n

        else:

            for i in range(n+1): t[i] = (1.0 - pow(f, i))/(1.0 - pow(f, n))

        if l: t *= l

        return t

    else:

        raise ValueError("n must be an integer, which is larger than 0")

# --------------------------------------------------------------------------- #
# --- End of module --------------------------------------------------------- #
# --------------------------------------------------------------------------- #

