#Main function for our program
from runge_kutta import RungeKutta as RK
from functions import *
import initialize
import matplotlib.pyplot as plt
import numpy as np

#Define all our constants. This is where we pass X, Y, Z
initialize.init(0.7,0.25,0.05)
val = np.logspace(-5., 5., 100)
count = 1
r_end = 10**10 #in units of m

for i in val:

    print i, "solar initial conditions"
    print count, "out of 100"
    count += 1
    # Initial core conditions of the star
    P0 = i*2.477 *10**16. # units of pascal
    M0 = 10**-4 # units of kg
    L0 = 10**-4 # units of J/s or W
    T0 = i*1.571*10**7 # units of K
    rho0 = i*1.622*10**5 #units of kg/m^3

    '''
    # Initial core conditions for the sun
    P0=2.477 *10**16. #in units of pascal
    M0=10**-4 #units of kg
    L0=10**-4 #units of J/s or W
    T0=1.571*10**7 #units of K
    rho0=1.622*10**5 #units of kg/m^3
    '''

    data = open('stars.dat','a')

    def stop(P,M,L,T,r,data):
	    if (P < data or T<data or r > 6.95*10**9.):
		    return True
	    else:
		    return False

    P,M,L,T,r,kappa_array,rho_array = RK(10**-6, r_end, 10.**5., P0, M0, L0, T0, rho0, dP, dM, dL, dT, rho, kappa, stop, 0.000001)
    #Need some way of recording the final luminosity and the final Temperature so we can put can plot these points onto an HR diagram

    #create file to write to 
    #write luminosity and temeprature to file each time we have a new set of parameters
    #print(len(r))
    '''
    plt.plot(r,T,label='T')
    plt.legend(loc='best')
    plt.show()
    plt.close()
     
    plt.plot(r,M,label='M')
    plt.legend(loc='best')
    plt.show()
    plt.close()
     
    plt.plot(r,L,label='L')
    plt.legend(loc='best')
    plt.show()
    plt.close()
     
    plt.plot(r,P,label='P')
    #plt.plot(r,L,label='L')
    #plt.xscale('log')
    #plt.yscale('log')
    #plt.ylim(0,10**25)
    plt.legend(loc='best')
    plt.show()

    '''

#    print(' optical depth:'+str(integrate(r,kappa_array,rho_array)))
#    print "radius", r[-1]
    print >> data, r[-1], M[-1], T[-1], L[-1]
    data.close()

 
