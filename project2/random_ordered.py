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
x=p
y=p
z=p
i=0
print(N)
while i<N:
	while x<(box_size-2*p):
		while y<(box_size-2*p):
			while z<(box_size-2*p):
				part.append(particle(m,[x,y,z],[0.,0.,0.]))
				z+=p*2
			y+=p*2
		x+=p*2
	N+=1

u=universe(N,0.,part)

u.writebinary('randomic_ordered.dat')

