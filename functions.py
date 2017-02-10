# I think I have done them incorrectly
# perhaps we even need another main function that includes this and runge_kutta to do the work

import scipy.constants
import numpy as np

X = 1
Y = 1
Z = 1
# WE'RE GOING TO NEED TO IMPORT X,Y,Z, AS WELL AS OTER INITIAL CONDITIONS
# define these globally

a = 7.565767 * 10 **(-16)
gamma = 5./3.
m_H = scipy.constants.m_p
k = scipy.constants.k
G = scipy.constants.G
c = scipy.constants.c
f_pp = 1.
f_3a = 1.
g_ff = 1.

mu = (2*X + 3./4. * Y + 1./2. * Z)**-1
# define non-changing constants
# these are defined globally


def T6(T):
   return T / 10**6

def T8(T):
   return T / 10**8
# define temperature constants


def rho(P,T):
   return (P - 1./3. * a * T**4) * mu * m_H / (k * T)
# define rho (and therefore mu)



def psi_pp(T):
   return 1 + 1.412 * 10 ** 8 * (1.0/X - 1) * np.exp(-49.98 * T6(T) ** (-1./3.))

def C_pp(T):
   return 1 + 0.0123 * T6(T)**(1./3.) + 0.0109 * T6(T)**(2./3.) + 0.000938 * T6(T)

def e_pp(P,T):
   return 0.241 * rho(P,T) * X**2 * f_pp * psi_pp(T) * C_pp(T) * T6(T)**(-2./3.) * np.exp(-33.8*T6(T)**(-1./3.))
# define epsilon_pp (and therefore psi_pp * C_pp)

def C_CNO(T):
   return 1 + 0.0027 * T6(T)**(1./3.) - 0.00778 * T6(T)**(2./3.) - 0.000149 * T6(T)
# define c_CNO

def e_CNO(P,T):
   return 8.67 * 10**20 * rho(P,T) * X * 0.5 * Z * C_CNO(T) * T6(T)**(-2./3.) * np.exp(-152.28*T6(T)**(-1./3.))
# define epsilon_CNO

def e_3a(P,T):
   return 50.9 * rho(P,T)**2 * Y**3 * T8(T)**-3 * f_3a * np.exp(-44.027 * T8(T)**-1)
# define epsilon_3alpha

def e(P,T):
   return e_pp(P,T) + e_CNO(T) + e_3a(P,T)
# define epsilon overall



def tdivgbf(P,T):
   return 0.708*(rho(P,T)*(1+X))**(1./5.)

def k_bf(P,T):
   return 4.34*10**21 * tdivgbf(P,T)**-1 * Z * (1+X) * rho(P,T) / T**3.5
# define kappa_bf (and therefore t/g_bf)

def k_ff(P,T):
   return 3.68*10**18 * g_ff * (1 - Z) * (1 + X) * rho(P,T) / T**3.5
# define kappa_ff

def k_es():
   return 0.02 * (1 + X)

def k_H(P,T):
   return 7.9*10**-34 * (Z / 0.02) * rho(P,T)**(1./2.) * T**9
# define kappa_H

def kappa(P,T):
   if (3000 <= T <= 6000) and (10**-7 <= rho(P,T) <= 10**-2) and (0.001 < Z < 0.003):
      return k_bf(P,T) + k_ff(P,T) + k_es + k_H(P,T)
   else:
      return k_bf(P,T) + k_ff(P,T) + k_es
# define kappa overall


# I've added ds so its clearer these are the differenetial eqautions

def dM(r,T,P,L):
   return 4 * np.pi * r**2 * rho(P,T)

def dP(r,T,M,L):
   return -G * M * rho(P,T) / r**2

def dL(r,T,M,P):
   return 4 * np.pi * r**2 * rho(P,T) * e(P,T)

# The left side of this equation needs to be a derivative
# But I'm not sure the best way to do so
def dT(r,P,M,L):
   if (log(P) / log(T) < (gamma / (gamma - 1.))):
      return -3. * kappa(P,T) * rho(P,T) * L / (4 * a * c * T**3 * 4 * np.pi * r**2)
   # radiative
   else:
      return -(1-1./gamma) * mu * m_H * G * M / (k * r**2)
   # adiabatic convection


