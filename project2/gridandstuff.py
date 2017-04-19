from nbody import *
from random import *
import numpy as np

#Read in initial conditions
box_size = 50
#u = universe.readbinary('randomic.dat')
u = universe.readbinary('init_L1.dat',scale=box_size)

#Declare parameters
Ncells = 2 * int(round((u.n)**(1/3),1))
omega_m = 0.3
omega_V = 0.7
omega_k = 1 - omega_m - omega_V
z = 50
a = 1.0 / (1+z)
da = 0.001
dx = box_size/Ncells
rho_0 = u.n/Ncells**3
count = 0

#Declare arrays
grid = np.zeros((Ncells,Ncells,Ncells))
phi_fft = np.zeros((Ncells,Ncells,Ncells),dtype=complex)

def G(l,m,n):
	#Function for calculating phi from delta
	if l == 0 and m == 0 and n == 0:
		return 0
	else:
		return -3 * omega_m/(8*a) * (1/(np.sin(np.pi*l/Ncells)**2+np.sin(np.pi*m/Ncells)**2+np.sin(np.pi*n/Ncells)**2))

def get_grid():
	#Calculates density of each cell
	grid = np.zeros((Ncells,Ncells,Ncells))
	for i in range(0,u.n):
		x,y,z = u.parts[i].position
		grid[int(x//dx)][int(y//dx)][int(z//dx)] += 1
	
def get_delta():
	#Calculate perturbations in density
	return (grid/rho_0) - 1

def get_phi(FFT):
	#Calculate the potential from delta
	for l in range(0,Ncells):
		for m in range(0,Ncells):
			for n in range(0,Ncells):	
				phi_fft[l][m][n] = FFT[l][m][n]*G(l,m,n)
	return np.fft.ifftn(phi_fft).real

def accels_grid():
	#Update acceleration of particles
	for N in range(0,u.n):
		u.parts[N].accel = [0.,0.,0.]
    
	for i in range(0,u.n):
		p_i = u.parts[i]
		x, y, z = p_i.position
		p_i.accel[0] = (-1.0/2.0) * (phi[(int(x//dx)+1) % Ncells][int(y//dx)][int(z//dx)] - phi[(int(x//dx)-1) % Ncells][int(y//dx)][int(z//dx)])
		p_i.accel[1] = (-1.0/2.0) * (phi[int(x//dx)][(int(y//dx)+1) % Ncells][int(z//dx)] - phi[int(x//dx)][(int(y//dx)-1) % Ncells][int(z//dx)])	
		p_i.accel[2] = (-1.0/2.0) * (phi[int(x//dx)][int(y//dx)][(int(z//dx)+1) % Ncells] - phi[int(x//dx)][int(y//dx)][(int(z//dx)-1) % Ncells])

def leapfrog_position_update(dt,mod=False,base=0):
	for i in range(0,u.n):
		p_i = u.parts[i]
		for k in range(3):
			p_i.position[k] += dt*p_i.velocity[k]/a**2
			if mod:
				p_i.position[k] = p_i.position[k] % base

def leapfrog_velocity_update(dt):
	for i in range(0,u.n):
		p_i = u.parts[i]
		for k in range(3):
			p_i.velocity[k] += 0.5*p_i.accel[k]*dt   

#Main loop
while a <= 1.0:
	#Calculate potential
	get_grid()
	phi = get_phi(np.fft.fftn(get_delta()))

	f = ((1/a)*(omega_m + omega_k*a + omega_V*a**3))**(-1/2)

	#Update position and velocity from potential
	dt=da*f
	leapfrog_position_update(dt,mod=True,base=box_size)
	leapfrog_velocity_update(dt)
	accels_grid()
	leapfrog_velocity_update(dt)
	u.write('./Data/universe{0:05d}.dat'.format(count))
	
	#Print current step and increment values
	print('Current a: ', round(a,3), end='\r')
	a += da
	count += 1
	
	





