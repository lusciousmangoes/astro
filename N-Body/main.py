from nbody import *

dt = 1./365.25
tf = 3.
t = 0.
u=universe.read("solar_system.dat")
u.accels()
print(u.T,u.U)

while t < tf:
	u.leapfrog_position_update(dt)
	u.leapfrog_velocity_update(dt)
	u.accels()
	u.leapfrog_velocity_update(dt)
	t += dt
	sun = u.parts[0]
	distances = ""
	for i in range(1,9):
		distances += str(u.parts[i].distance(sun))+ " " 
	print(t,u.T+u.U,distances)

