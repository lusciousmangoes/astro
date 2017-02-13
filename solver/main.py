#Main function for our program
from runge_kutta import RungeKutta as RK
from functions import *

#Taking zero to be the center of the star
#initial pressure at center of star
P0=1.0
#initial mass at center of star
M0=1.0
#initial luminosity at center of star
L0=1.0
#initial temperature at center of star
T0=1.0
#initial density at the center of the star
rho0=1.0

r_end

def stop(P,M,L,T,data):
	if (P < data or M<data or L<data or T<data):
		return True
	else:
		return False		

RK(0.0, r_end, 0.01, P0, M0, L0, T0, rho0, dP, dM, dL, dT, stop, 0.00001):

 
