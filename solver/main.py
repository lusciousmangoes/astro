#Main function for our program
from runge_kutta import RungeKutta as RK
from functions import *
import initialize
import matplotlib.pyplot as plt

#Define all our constants. This is where we pass X, Y, Z, and it should work everywhere
initialize.init(0.7,0.25,0.05)
data = open('stars.dat','a')

#Taking zero to be the center of the star
#initial pressure at center of star
P0=2.477 *10**17. #in units of pascal
#P0=1.0
#initial mass at center of star
M0=10**-4 #units of kg

#initial luminosity at center of star
L0=10**-4 #units of J/s or W
#L0=1.0
#initial temperature at center of star
T0=1.571*10**8 #units of K
#T0=1.0
#initial density at the center of the star
rho0=1.622*10**5 #units of kg/m^3
#rho0=1.0
#final r value, file should end though with stop function
r_end = 10**10 #in units of m
#r_end=10.0
def stop(P,M,L,T,r,data):
	if (P < data or T<data or r > 6.95*10**9.):
                print(T,P,r)
		return True
	else:
		return False

#print e(10.0,10.0**7) #Testing some stuff	

#print(dP(0.01,T0,M0,P0))
P,M,L,T,r,kappa_array,rho_array = RK(10**-6, r_end, 10.**5., P0, M0, L0, T0, rho0, dP, dM, dL, dT, rho, kappa, stop, 0.000001)
#Need some way of recording the final luminosity and the final Temperature so we can put can plot these points onto an HR diagram

#create file to write to 
#write luminosity and temeprature to file each time we have a new set of parameters
#print(len(r))
plt.plot(r,T,label='T')
plt.legend(loc='best')
#plt.show()
 
plt.plot(r,M,label='M')
plt.legend(loc='best')
#plt.show()
 
plt.plot(r,L,label='L')
plt.legend(loc='best')
#plt.show()
 
plt.plot(r,P,label='P')
#plt.plot(r,L,label='L')
#plt.xscale('log')
#plt.yscale('log')
#plt.ylim(0,10**25)
plt.legend(loc='best')
#plt.show()

print(' optical depth:'+str(integrate(r,kappa_array,rho_array)))

print >> data, r[-1], M[-1], T[-100], L[-1]
data.close()




 
