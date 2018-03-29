#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# October 2016
# Pascal Beckstein (p.beckstein@hzdr.de)

L = 4.00 * 1e-3             # [L]        = m
R = 0.15 * 1e-3             # [R]        = m
nu = 1e-6                   # [nu]       = m^2/s
rho = 998.00                # [rho]      = kg/m^3
eta = nu*rho                # [eta]      = m^2/s * kg/m^3 = kg/s/m
#ddxT = 25.4e+3              # [ddxT]     = K/m
ddxT = 25.4e+3/500          # [ddxT]     = K/m
ddTSigma = 0.0274e-3        # [ddTSigma] = N/m/K
alpha = 1.08e-9             # [alpha]    = m^2/s


u_0 = R*ddxT/eta * ddTSigma # [u_0]      = m*K/m * s*m/kg * kg*m/m/K/s^2 = m/s
u_inf = 0.5 * u_0

Re = u_inf*R / nu
Mg = ddTSigma * R*(R*ddxT) / (eta*alpha)

print("R            = %g m" % (R))
print("nu           = %g m^2/s" % (nu))
print("rho          = %g kg/m^3" % (rho))
print("eta          = %g kg/s/m" % (eta))
print("ddxT         = %g K/m" % (ddxT))
print("ddTSigma     = %g N/m/K" % (ddTSigma))


print("u_0          = %g m/s = %g mm/s" % (u_0,u_0*1e+3))
print("u_inf        = %g m/s = %g mm/s" % (u_inf,u_inf*1e+3))

print("Re           = %g" % (Re))
print("Mg           = %g" % (Mg))