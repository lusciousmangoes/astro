import scipy.constants
import numpy as np

X = 1
Y = 1
Z = 1
# WE'RE GOING TO NEED TO IMPORT X,Y,Z, AS WELL AS OTER INITIAL CONDITIONS


a = 7.565767 * 10 **(-16)
m_H = scipy.constants.m_p
k = scipy.constants.k
f_pp = 1.
f_3a = 1.


def T6(T):
   val = T / 10**6
   return val

def T8(T):
   val = T / 10**8
   return val

def C_pp(T):
   val = 1 + 0.0123 * T6(T)**(1./3.) + 0.0109 * T6(T)**(2./3.) + 0.000938 * T6(T)
  return val

def mu (X,Y,Z):
   val = (2*X + 3./4. * Y + 1./2. * Z)**-1
   return val

def rho(P,T):
   val = (P - 1./3. * a * T**4) * mu(X,Y,Z) * m_H / (k * T)
   return val

def psi_pp(X,T):
   val = 1 + 1.412 * 10 ** 8 * (1/X - 1) * np.exp(-49.98 * T6(T) ** (-1./3.))
   return val

def e_pp(P,T):
   val = 0.241 * rho(P,T) * X**2 * f_pp * psi_pp(X,T) * C_pp(T) * T6(T)**(-2./3.) * np.exp(-33.8*T6(T)**(-1./3.))
   return val

def e_CNO(T):
   val = 1 + 0.0027 * T6(T)**(1./3.) - 0.00778 * T6(T)**(2./3.) - 0.000149 * T6(T)
   return val

def e_3a(
   val = 50.9 * rho(P,T)**2 * Y**3 * T8(T)**-3 * f_3a * np.exp(-44.027 * T8(T)**-1)
   return val




'''
m

epsilon

P

L

T

kappa
'''



