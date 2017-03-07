import numpy as np
import initialize as cons

def T6(T):
   return T / 10.0**6

def T8(T):
   return T / 10.0**8
# define temperature constants

def rho(P,T):
   return abs((P - 1./3. * cons.a * T**4) * cons.mu * cons.m_H / (cons.k * T))
# define rho

def psi_pp(T):
   return 1 + ((1.412 * 10.0 ** 8.0) * ((1.0/cons.X) - 1)) * np.exp(-49.98 * T6(T) ** (-1./3.))

def C_pp(T):
   return 1 + 0.0123 * T6(T)**(1./3.) + 0.0109 * T6(T)**(2./3.) + 0.000938 * T6(T)

def e_pp(P,T):
   return 0.241 * rho(P,T) * cons.X**2 * cons.f_pp * psi_pp(T) * C_pp(T) * T6(T)**(-2./3.) * np.exp(-33.8*T6(T)**(-1./3.))
# define epsilon_pp (and therefore psi_pp * C_pp)

def C_CNO(T):
   return 1 + 0.0027 * T6(T)**(1./3.) - 0.00778 * T6(T)**(2./3.) - 0.000149 * T6(T)
# define c_CNO

def e_CNO(P,T):
   return 8.67 * 10**20 * rho(P,T) * cons.X * 0.5 * cons.Z * C_CNO(T) * T6(T)**(-2./3.) * np.exp(-152.28*T6(T)**(-1./3.))
# define epsilon_CNO

def e_3a(P,T):
   return 50.9 * rho(P,T)**2 * cons.Y**3 * T8(T)**-3 * cons.f_3a * np.exp(-44.027 * T8(T)**-1)
# define epsilon_3alpha

def e(P,T):
   return e_pp(P,T) + e_CNO(P,T) + e_3a(P,T)
# define epsilon overall



def tdivgbf(P,T):
   return 0.708*(abs(rho(P,T))*(1+cons.X))**(1./5.)

def k_bf(P,T):
   return 4.34*10**21 * tdivgbf(P,T)**-1 * cons.Z * (1+cons.X) * rho(P,T) / T**3.5
# define kappa_bf (and therefore t/g_bf)

def k_ff(P,T):
   return 3.68*10**18 * cons.g_ff * (1 - cons.Z) * (1 + cons.X) * rho(P,T) / T**3.5
# define kappa_ff

def k_es():
   return 0.02 * (1 + cons.X)

def k_H(P,T):
   return 7.9*10**-34 * (cons.Z / 0.02) * rho(P,T)**(1./2.) * T**9
# define kappa_H

def kappa(P,T):
   if (3000 <= T <= 6000) and (10**-7 <= rho(P,T) <= 10**-2) and (0.001 < cons.Z < 0.003):
      return k_bf(P,T) + k_ff(P,T) + k_es() + k_H(P,T)
   else:
      return k_bf(P,T) + k_ff(P,T) + k_es()
# define kappa overall


# I've added ds so its clearer these are the differenetial eqautions

def dM(r,T,P,L):
   return 4 * np.pi * r**2 * rho(P,T)

def dP(r,T,M,P):
   return -cons.G * M * rho(P,T) / r**2

def dL(r,T,M,P):
   return 4 * np.pi * r**2 * rho(P,T) * e(P,T)

def dT(r,P,M,L,T):
   if (16 * cons.G * M * cons.a * cons.c * T**4 * np.pi / (3 * P * kappa(P,T) * L) < (cons.gamma / (cons.gamma - 1.))):
      #print "convective"
      return -(1-1./cons.gamma) * cons.mu * cons.m_H * cons.G * M / (cons.k * r**2)
   # adiabatic convection
   else:
      #print "radiative"
      return -3. * kappa(P,T) * rho(P,T) * L / (4 * cons.a * cons.c * T**3 * 4 * np.pi * r**2)
   # radiative  
    
def integrate(r,kappa,rho):
   A = 0.0 
   integrand=kappa*rho
   for i in range(0,len(r)-1):
      A += ((integrand[i]+integrand[i+1])/2.0)*(r[i+1]-r[i])
   return A

# I've fixed the "if" conditions, so now main runs without any error
