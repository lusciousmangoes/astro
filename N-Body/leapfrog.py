import numpy as np

def integrator():
	accels()

	for i in range N
		# Compute values at next time step
		vh = v[0][i] + 0.5*a[0][i]*dt;
		r[0][i] = r[0][i] + vh*dt;
		v[0][i] = vh + 0.5*a[0][i]*dt;
	
		vh = v[1][i] + 0.5*a[1][i]*dt;
		r[1][i] = r[1][i] + vh*dt;
		v[1][i] = vh + 0.5*a[1][i]*dt;
	
		vh = v[2][i] + 0.5*a[2][i]*dt;
		r[2][i] = r[2][i] + vh*dt;
		v[2][i] = vh + 0.5*a[2][i]*dt;

