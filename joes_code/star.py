from __future__ import print_function
import numpy as np
from physics import *
import sys

# take 1.  Set P = PO, guess at T, and integrate outward.

X = 0.73
Y = 0.25
Z = 1.0 - X - Y


def dPdr(r, P, M, L, T):
    rho = density(P, T, X, Y, Z)
    return -G * M * rho / r / r

def dMdr(r, P, M, L, T):
    rho = density(P, T, X, Y, Z)
    return 4.0 * np.pi * r * r * rho
     
def dLdr(r, P, M, L, T):
    rho = density(P, T, X, Y, Z)
    e = epsilon(rho, T, X, Y, Z)    
    return 4.0 * np.pi * r * r * rho * e
    
def dTdr_rad(r, P, M, L, T):
    rho = density(P, T, X, Y, Z)
    kappa = opacity(rho, T, X, Y, Z)    
    return -3.0 / 4.0 / a_rad / c * kappa * rho / T / T / T * L / 4.0 / np.pi / r / r
    
def dTdr_con(r, P, M, L, T):
    mu = molecular_weight(X, Y, Z)
    return -1.0 * (1.0 - 1.0/gamma) * mu * m_H / k_B * G * M / r / r


dr = 1e4
P0 = 2.342e16
T0 = 1.57e7

r = np.arange(1e4, 1e9, dr)
N = len(r)
P = np.zeros(N); P[0] = P0
M = np.zeros(N)
L = np.zeros(N)
T = np.zeros(N); T[0] = T0

dT = 0.0
dP = 0.0
ratio = 3.0
transport_type = "r"
for i in range(N-1):
    kP1 = dr * dPdr(r[i], P[i], M[i], L[i], T[i])
    kM1 = dr * dMdr(r[i], P[i], M[i], L[i], T[i])
    kL1 = dr * dLdr(r[i], P[i], M[i], L[i], T[i])
    
    dTdr_r = dTdr_rad(r[i], P[i], M[i], L[i], T[i])
    dTdr_c = dTdr_con(r[i], P[i], M[i], L[i], T[i])
    
    if i > 1:
        dT = T[i] - T[i-1]
        dP = P[i] - P[i-1]
        ratio = (T[i] + T[i-1]) / (P[i] + P[i-1]) * dP / dT
        if ratio < gamma / (gamma - 1.0):
            kT1 = dr * dTdr_c
            transport_type = "c"
        else:
            kT1 = dr * dTdr_r
            transport_type = "r"
    else:
        kT1 = dr * dTdr_r
    
    #forget it, just do a fully convective star, it look nicer
    #kT1 = dr * dTdr_c 
            
    P[i+1] = P[i] + kP1
    M[i+1] = M[i] + kM1
    L[i+1] = L[i] + kL1
    T[i+1] = T[i] + kT1
    
    if T[i+1] < 0.0 or P[i+1] < 0.0:
        print("Terminating integration at r = %5.3f (T = %5.3g, P = %5.3g)" % (r[i+1]/R_sun,  T[i+1], P[i+1]), file=sys.stderr)
        break
     
    rho = density(P[i+1], T[i+1], X, Y, Z)
    kappa = opacity(rho, T[i+1], X, Y, Z)
    
    print(i+1, r[i+1]/R_sun, P[i+1]/1e16, M[i+1]/M_sun, L[i+1]/L_sun, T[i+1], rho, ratio, kappa)
    print("Zone %0d:  transport is %s (dTdr_r = %5.3f, dTdr_c = %5.3f, ratio=%5.3f)" % (i+1, transport_type, dTdr_r, dTdr_c, ratio), file=sys.stderr)

