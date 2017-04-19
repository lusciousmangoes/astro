from nbody import *
from random import *
import numpy as np

box_size = 50
Ncells=16
N=Ncells**3
part=[]
H0 = 73.8 *((3.086**-1)*10**-19)
G = 6.674*10**-11
rho = 0.27*(3*H0**2)/(8*np.pi*G)
V=(box_size)*(box_size)*(box_size)
m = rho*V/N
p = box_size/Ncells/2


x=np.linspace(0,box_size,Ncells)
y=np.linspace(0,box_size,Ncells)
z=np.linspace(0,box_size,Ncells)

print x


for i in range(len(x)):
	for j in range(len(y)):
		for k in range(len(z)):
			part.append(particle(m,[x[i],y[j],z[k]],[0.,0.,0.]))



u=universe(N,0.,part)

u.write('randomic_ordered.dat')

