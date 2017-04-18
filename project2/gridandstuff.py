from nbody import *
from random import *
import numpy as np

a=1/(1+50)
box_size = 50
u=universe.readbinary('randomic.dat')
Ncells=16
grid=np.zeros((Ncells,Ncells,Ncells))
delta=np.zeros((Ncells,Ncells,Ncells))
phi_fft=np.zeros((Ncells,Ncells,Ncells),dtype=complex)
gx=np.zeros((Ncells,Ncells,Ncells))
gy=np.zeros((Ncells,Ncells,Ncells))
gz=np.zeros((Ncells,Ncells,Ncells))
dx=box_size/Ncells
rho_0 = u.n/Ncells**3
omega = 0.27

for i in range(0,u.n):
	x,y,z = u.parts[i].position
	grid[int(x//dx)][int(y//dx)][int(z//dx)]+=1
	
#print(sum(sum(sum(grid))))
delta=(grid/rho_0) - 1
print(sum(sum(sum(delta))))

delta_fft=np.fft.fftn(delta)
def G(l,m,n):
	if l==0 and m==0 and n==0:
		return 0
	else:
		return -3*omega/(8*a) * (1/(np.sin(np.pi*l/Ncells)**2+np.sin(np.pi*m/Ncells)**2+np.sin(np.pi*n/Ncells)**2))

for l in range(0,Ncells):
		for m in range(0,Ncells):
			for n in range(0,Ncells):	
				phi_fft[l][m][n] = delta_fft[l][m][n]*G(l,m,n)
phi=np.fft.ifftn(phi_fft)
phi = phi.real

print(phi.shape)

def accels_grid():
	for N in range(0,u.n):
		u.parts[N].accel=[0.,0.,0.]
    
	for i in range(0,u.n):
		p_i = u.parts[i]
		x, y, z = p_i.position
		#print((int(x//dx)+1) % Ncells)
		p_i.accel[0] = (-1.0/2.0) * (phi[(int(x//dx)+1) % Ncells][int(y//dx)][int(z//dx)] - phi[(int(x//dx)-1) % Ncells][int(y//dx)][int(z//dx)])
		p_i.accel[1] = (-1.0/2.0) * (phi[int(x//dx)][(int(y//dx)+1) % Ncells][int(z//dx)] - phi[int(x//dx)][(int(y//dx)-1) % Ncells][int(z//dx)])	
		p_i.accel[2] = (-1.0/2.0) * (phi[int(x//dx)][int(y//dx)][(int(z//dx)+1) % Ncells] - phi[int(x//dx)][int(y//dx)][(int(z//dx)-1) % Ncells])

while z>1:
	a=1/(1+z)
	da=1
	u.leapfrog_position_update(da)
	u.leapfrog_velocity_update(da)
	accels_grid()
	u.leapfrog_velocity_update(da)
	z-=1
	
	





