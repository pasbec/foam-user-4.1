#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Compute the discrete fourier transform of data from file
# Jan 2018
# Vladimir Galindo (v.galindo@hzdr.de)
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

import numpy as np
import pylab as p

# --------------------------------------------------------------------------- #
# --- Parameters ------------------------------------------------------------ #
# --------------------------------------------------------------------------- #

print(sys.argv)

if len(sys.argv) > 1:
    filename = sys.argv[1]
else:
    print("Missing file name.")
    sys.exit(1)

if len(sys.argv) > 2:
    deltaT = float(sys.argv[2])
else:
    print("Missing time step width.")
    sys.exit(1)

# --------------------------------------------------------------------------- #
# --- Main program sequence ------------------------------------------------- #
# --------------------------------------------------------------------------- #

with open(filename) as file:
     floats = map(float,file)

# data
UX = np.array(floats)
# number of data values
NT = len(floats)

print("data average       : %f" % (UX.mean()))
UX = UX - UX.mean()

# time step DT = (tmax-tmin)/(NT-1)
DT = deltaT
print("delta_t            : ",DT)

# smallest frequency
DF = 1.0/((NT-1.0)*DT)
print("smallest frequency : %f" % (DF))

# real fourier transform
fft=(2.0/NT)*np.fft.rfft(UX)
# power spectrum
power=abs(fft)
omit=3
maxPowerID=np.argmax(power[omit:])+omit
if NT%2==0:
   F = np.linspace(0,NT*DF/2,NT/2+1)
else:
   F = np.linspace(0,(NT-1)*DF/2,(NT-1)/2+1)

print("N                  : %d" % (NT))
print("peak frequency     : %f Hz" % (F[maxPowerID]))
print("Periode            : %f s"  % (1.0/F[maxPowerID]))

p.plot(F[omit:NT/2],power[omit:NT/2])
p.axvline(F[maxPowerID],c="r")
p.annotate(" f = %g Hz, T = %g s" % (F[maxPowerID], (1.0/F[maxPowerID])), (F[maxPowerID]*1.03,power[maxPowerID]*0.95))
p.xlabel("$f / Hz$")
p.ylabel("power spectrum: |rfft(data)|")
p.xlim(0,F[NT/2/5])
p.grid()

p.savefig(os.path.dirname(filename) + "/" + os.path.splitext(os.path.basename(filename))[0] + ".pdf")

# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
# --------------------------------------------------------------------------- #
