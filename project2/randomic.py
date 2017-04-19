from nbody import *
from random import *
import numpy as np

box_size = 50
N=32**3
H0 = 73.8 *((3.086**-1)*10**-19)
G = 6.674*10**-11
rho = 0.27*(3*H0**2)/(8*np.pi*G)
V=50*50*50
m = rho*V/N
part=[]
i=0
print(N)
while i<N :
    a = uniform(0,box_size)
    b = uniform(0,box_size)
    c = uniform(0,box_size)
    part.append(particle(m,[a,b,c],[0.,0.,0.]))
    i+=1

u=universe(N,0.,part)

u.writebinary('randomic.dat')

