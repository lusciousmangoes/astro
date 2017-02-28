import numpy as np
from nbody import particle
G=1
earth = particle(0.1,[0.9,0,0],[0,-1.2,0])
sun = particle(0.9,[-0.1,0,0],[0,0.133,0])

E0 = 0.5 *earth.mass * earth.vel2() + 0.5 * sun.mass*sun.vel2() -G *sun.mass*earth.mass/earth.distance(sun) 


dt = 0.001
for t in np.arange(0,600,dt):
   r3 = earth.distance(sun)**3.
   
   for i in range(3):
      earth.accel[i]= -G*sun.mass*(earth.position[i]-sun.position[i])/r3
      sun.accel[i] = -G*earth.mass*(sun.position[i]-earth.position[i])/r3
      earth.velocity[i] += earth.accel[i]*dt 
      sun.velocity[i] += sun.accel[i]*dt

      earth.position[i] += earth.velocity[i]*dt
      sun.position[i] += sun.velocity[i]*dt
   E = 0.5 *earth.mass * earth.vel2() + 0.5 * sun.mass*sun.vel2() -G *sun.mass*earth.mass/earth.distance(sun) 
   if(t/dt) % 1000 == 0:
      print t,earth.x(),earth.y(),earth.z(),sun.x(),sun.y(),sun.z(),(E-E0)/E0 
