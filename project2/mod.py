import numpy as np
from nbody import *


u = universe.readbinary('init_L1.dat')

H=14722.314
#Mpc/km/s




x=np.zeros(u.n)
y=np.zeros(u.n)
z=np.zeros(u.n)
vx=np.zeros(u.n)
vy=np.zeros(u.n)
vz=np.zeros(u.n)
flag=np.zeros(u.n)


for i in range(0,u.n):
	x[i]=50*u.parts[i].position[0]-25
	y[i]=50*u.parts[i].position[1]-25
	z[i]=50*u.parts[i].position[2]-25
#	vx[i]=H*u.parts[i].position[0]
#	vy[i]=H*u.parts[i].position[1]
#	vz[i]=H*u.parts[i].position[2]
	vx[i]=H*x[i]
	vy[i]=H*y[i]
	vz[i]=H*z[i]




for i in range(u.n):
	if (x[i]**2+y[i]**2+z[i]**2)>625:
		flag[i]=1


numba=int(u.n-sum(flag))
masss=1./u.n




name="./DNC/changed.dat"

f = open(name,"wb")
f.write(pack('<i',numba))
f.write(pack('<d',u.t))
a=array.array('f')
for i in range(0,u.n):
	if flag[i] == 0:
		a.fromlist([masss,x[i],y[i],z[i],u.parts[i].velocity[0],u.parts[i].velocity[1],u.parts[i].velocity[2]])
a.tofile(f)

f.close()






