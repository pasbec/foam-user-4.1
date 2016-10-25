#!/usr/bin/python
# October 2016
# Vladimir Galindo (v.galindo@hzdr.de)
# Pascal Beckstein (p.beckstein@hzdr.de)

# --------------------------------------------------------------------------- #
# --- Libraries ------------------------------------------------------------- #
# --------------------------------------------------------------------------- #

import math
import numpy as np
import scipy as sp
from scipy import special as spsp

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

N = 40
H0 = 0.03
R0 = 0.03
sigma = 3.2906e+06
omega = 2.0*math.pi*50.0
B0 = 0.421599e-3

nr = 61
nz = 61

print "sigma      : %g S/m" % (sigma)
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
            a = a + (l[i]*math.cosh(l[i]*z))/((l[i]*l[i]-1.0)*spsp.jv(1,l[i])*math.cosh(l[i]*H0/R0))
        b = r - a*r
    else:
        a = 0.0
        for i in range(len(l)):
            a = a + (spsp.jv(1,l[i]*r)*math.cosh(l[i]*z))/((l[i]*l[i]-1.0)*spsp.jv(1,l[i])*math.cosh(l[i]*H0/R0))
        b = r - 2.0*a
    return b



# Lorentz-Force density
def fl(r,z):

    return 0.5*sigma*omega*B0*B0*R0*f(r/R0,z/R0)

# --------------------------------------------------------------------------- #
# --- Main program sequence ------------------------------------------------- #
# --------------------------------------------------------------------------- #

# Mesh and discretisation
rl = np.linspace( 0.0,R0,nr)
zl = np.linspace(0.0,H0,nz)

R,Z = np.meshgrid(rl,zl,indexing='ij')
F   = np.zeros(R.shape)


# Calculate and write data

fileName = 'LorentzForce.dat'

with open(fileName,'w') as lff:

    lff.write('# Variables: ' + 'r' + ' ' + 'z' + ' ' + 'F' + ' \n')

    lff.write('# Resolution: ' + 'nr=' + str(nr) + ', ' + 'nz=' + str(nz) + ' \n')

    for i in range(len(rl)):
        for j in range(len(zl)):

            F[i,j] = fl(rl[i],zl[j])

            lff.write(str(rl[i]) + ' ' + str(zl[j]) + ' ' + str(F[i,j]) + ' \n')

print "Fmin, Fmax :", str(np.min(F)) + ",", str(np.max(F))

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
