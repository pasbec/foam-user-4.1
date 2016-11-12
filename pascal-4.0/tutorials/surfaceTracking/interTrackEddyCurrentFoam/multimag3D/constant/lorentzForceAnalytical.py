#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
# Vladimir Galindo (v.galindo@hzdr.de)
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import math as m
import numpy as np
from scipy import special as spsp

from foamTools.ioInfo import fileGetPath, objectIndent, objectHeader, objectFooter

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

N = 40
H0 = 0.03
R0 = 0.03
sigma = 3.2906e+06
omega = 2.0*m.pi*50.0
B0 = 0.421599e-3

nr = 61
nz = 61

print "N          : %g" % (N)
print "H0, R0     :", str(H0) + ",", str(R0)
print "sigma      : %g S/m" % (sigma)
print "omega      : %g" % (omega)
print "B0         : %g T" % (B0)

print "nr, nz     :", str(nr) + ",", str(nz)

# --------------------------------------------------------------------------- #
# --- Function definitions -------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Roots of Bessel-fct
l = spsp.jnp_zeros(1,N)



# Base functions for RMF
def f(r,z):

    if (r==0):
        a = 0.0
        for i in range(len(l)):
            a = a + (l[i]*m.cosh(l[i]*z))/((l[i]*l[i]-1.0)*spsp.jv(1,l[i])*m.cosh(l[i]*H0/R0))
        b = r - a*r
    else:
        a = 0.0
        for i in range(len(l)):
            a = a + (spsp.jv(1,l[i]*r)*m.cosh(l[i]*z))/((l[i]*l[i]-1.0)*spsp.jv(1,l[i])*m.cosh(l[i]*H0/R0))
        b = r - 2.0*a
    return b



# Lorentz-Force density
def fl(r,z):

    return 0.5*sigma*omega*B0*B0*R0*f(r/R0,z/R0)

# --------------------------------------------------------------------------- #
# --- Main program sequence ------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Mesh and discretisation
rl = np.linspace(0.0,R0,nr)
zl = np.linspace(0.0,H0,nz)

R,Z = np.meshgrid(rl,zl,indexing="ij")
F   = np.zeros(R.shape)


# Calculate and write data
with open(fileGetPath("lorentzForceAnalytical.dat"),"w") as lff:

    lff.write("# Variables: " + "r" + " " + "z" + " " + "F" + "\n")

    lff.write("# Resolution: " + "nr=" + str(nr) + ", " + "nz=" + str(nz) + "\n")

    for ri in range(nr):
        for zi in range(nz):

            F[ri,zi] = fl(rl[ri],zl[zi])

            lff.write(str(rl[ri]) + " " + str(zl[zi]) + " " + str(F[ri,zi]) + "\n")

print "Fmin, Fmax :", str(np.min(F)) + ",", str(np.max(F))

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
