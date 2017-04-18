from nbody import *
from random import *
import numpy as np

#Read in initial conditions
u = universe.readbinary('randomic.dat')

#Declare parameters
box_size = 50
Ncells = 16
omega = 0.27
z = 50
a = 1.0 / (1+z)
da = 0.01
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
		return -3 * omega/(8*a) * (1/(np.sin(np.pi*l/Ncells)**2+np.sin(np.pi*m/Ncells)**2+np.sin(np.pi*n/Ncells)**2))

def get_grid():
	#Calculates density of each cell
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

#Main loop
while a <= 1.0:
	get_grid()
	phi = get_phi(np.fft.fftn(get_delta())) #Calculate phi from the Fourier transform of delta

	u.leapfrog_position_update(da)
	u.leapfrog_velocity_update(da)
	accels_grid()
	u.leapfrog_velocity_update(da)
	u.write('./Data/universe{0:05d}.dat'.format(count))
	
	print('Current time: ', a)
	a += da
	count += 1
	
	





