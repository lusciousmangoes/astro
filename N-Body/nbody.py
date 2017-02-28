from math import sqrt
class particle:
   def __init__(self, mass, position,velocity):
      self.mass = mass
      self.position = position
      self.velocity = velocity
      self.accel =[0,0,0]

   def x(self):
      return self.position[0]
   def y(self):
      return self.position[1] 
   def z(self):
      return self.position[2]


   def distance(self,other):
      return sqrt((self.x()-other.x())**2+(self.y()-other.y())**2+(self.z()-other.z())**2)   

   def vel2(self):
      return self.velocity[0]**2. + self.velocity[1]**2. + self.velocity[2]**2.

