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
# define non-changing constants
# these are defined globally

print G

def T6(T):
   val = T / 10**6
   return val

def T8(T):
   val = T / 10**8
   return val
# define temperature constants


def mu (X,Y,Z):
   val = (2*X + 3./4. * Y + 1./2. * Z)**-1
   return val

def rho(P,T):
   val = (P - 1./3. * a * T**4) * mu(X,Y,Z) * m_H / (k * T)
   return val
# define rho (and therefore mu)



def psi_pp(T):
   val = 1 + 1.412 * 10 ** 8 * (1/X - 1) * np.exp(-49.98 * T6(T) ** (-1./3.))
   return val

def C_pp(T):
   val = 1 + 0.0123 * T6(T)**(1./3.) + 0.0109 * T6(T)**(2./3.) + 0.000938 * T6(T)
   return val

def e_pp(P,T):
   val = 0.241 * rho(P,T) * X**2 * f_pp * psi_pp(T) * C_pp(T) * T6(T)**(-2./3.) * np.exp(-33.8*T6(T)**(-1./3.))
   return val
# define epsilon_pp (and therefore psi_pp * C_pp)

def e_CNO(T):
   val = 1 + 0.0027 * T6(T)**(1./3.) - 0.00778 * T6(T)**(2./3.) - 0.000149 * T6(T)
   return val
# define epsilon_CNO

def e_3a(P,T):
   val = 50.9 * rho(P,T)**2 * Y**3 * T8(T)**-3 * f_3a * np.exp(-44.027 * T8(T)**-1)
   return val
# define epsilon_3alpha

def e(P,T):
   val = e_pp(P,T) + e_CNO(T) + e_3a(P,T)
   return val
# define epsilon overall



def tdivgbf(P,T):
   val = 0.708*(rho(P,T)*(1+X))**(1./5.)

def k_bf(P,T):
   val = 4.34*10**21 * tdivgbf(P,T)**-1 * Z * (1+X) * rho(P,T) / T**(3./5.)
   return val
# define kappa_bf (and therefore t/g_bf)

def k_ff(P,T):
   val = 3.68*10**18 * g_ff * (1 - Z) * (1 + X) * rho(P,T) / T**(3./5.)
   return val
# define kappa_ff

def k_es():
   val = 0.02 * (1 + X)
   return val
# define kappa_es

def k_H(P,T):
   val = 7.9*10**-34 * (Z / 0.02) * rho(P,T)**(1./2.) * T**9
# define kappa_H

def kappa(P,T):
   if (3000 <= T <= 6000) and (10**-7 <= rho(P,T) <= 10**-2) and (0.001 < Z < 0.003):
      val = k_bf(P,T) + k_ff(P,T) + k_es + k_H(P,T)
   else:
      val = k_bf(P,T) + k_ff(P,T) + k_es
   return val
# define kappa overall






# below here i'm pretty sure there are errors
# quite possibly in the if statements as well

def M(r,T,P,L):
   val = 4 * np.pi * r**2 * rho(P,T)
   return val

def P(r,T,M,L):
   val = -G * M rho(P,T) / r**2
   return val

def L(r,T,M,P):
   val = 4 * np.pi * r**2 * rho(P,T) * e(P,T)
   return val

def T(r,P,M,L):
   if (log(P) / log(T) < (gamma / (gamma - 1)):
      val = -3. * kappa(P,T) * rho(P,T) * L / (4 * a * c * T**3 * 4 * np.pi * r**2)
   # radiative
   else:
      val = -(1-1/gamma) * mu(X,Y,Z) * m_H * G * M / (k * r**2)
   # adiabatic convection


