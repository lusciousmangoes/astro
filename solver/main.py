#Main function for our program
from runge_kutta import RungeKutta as RK
from functions import *
import initialize
import matplotlib.pyplot as plt

#Define all our constants. This is where we pass X, Y, Z, and it should work everywhere
initialize.init(0.7,0.25,0.05)

#Taking zero to be the center of the star
#initial pressure at center of star
P0=247700 #in units of pascal
#P0=1.0
#initial mass at center of star
M0=10**-4 #units of kg

#initial luminosity at center of star
L0=10**-4 #units of J/s or W
#L0=1.0
#initial temperature at center of star
T0=1.571*10**7 #units of K
#T0=1.0
#initial density at the center of the star
rho0=1.622*10**5 #units of kg/m^3
#rho0=1.0
#final r value, file should end though with stop function
r_end = 10**2 #in units of m
#r_end=10.0
def stop(P,M,L,T,data):
	if (P < data or M<data or L<data or T<data):
		return True
	else:
		return False

#print e(10.0,10.0**7) #Testing some stuff	

P,M,L,T,r = RK(0.01, r_end, 0.01, P0, M0, L0, T0, rho0, dP, dM, dL, dT, stop, 0.0001)

print(L)
#Need some way of recording the final luminosity and the final Temperature so we can put can plot these points onto an HR diagram

#create file to write to 
#write luminosity and temeprature to file each time we have a new set of parameters

plt.plot(r,M,label='M')
plt.plot(r,L,label='L')
plt.legend(loc='best')
plt.show()

 
