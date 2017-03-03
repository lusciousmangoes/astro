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


   def __str__(self):
      return "{0:.10e} {1:.10e} {2:.10e} {3:.10e} {4:.10e} {5:.10e} {6:.10e}".format(self.mass,self.position[0],self.position[1],self.position[2],self.velocity[0],self.velocity[1],self.velocity[2])


   @classmethod
   def string_thing(cls,s):
      var = s.split()
      mass = float(var[0])
      x = float(var[1])
      y = float(var[2])
      z = float(var[3])
      vx = float(var[4])
      vy = float(var[5])
      vz = float(var[6])
      return cls(mass,[x,y,z],[vx,vy,vz])


class universe:
   def __init__(self, n, t, parts, G=1.):
      self.n = n
      self.t = t
      self.parts = parts
      self.G = G
      self.T = 0.
      self.U = 0.



#   def accels():

   @classmethod
   def read(cls,fname):
      f = open(fname, "r")
      n = int(f.readline())
      t = float(f.readline())
      G = float(f.readline())
      parts = []
      for i in range(0,n):
         parts.append(particle.string_thing(f.readline()))
      return cls(n,t,parts,G)



   def write(self, name):
      f = open(name,"w")
      f.write("{0}\n{1}\n{2}\n".format(self.n,self.t,self.G))
      for i in range(0,self.n):
         f.write(str(self.parts[i])+"\n")
      f.close()




def main():
#   part = particle(1,[0.2,0.2,0.2],[6.,6.,6.])
#   u = universe(1,5,[part],7.)
#   u.write("fileIO.txt")

   u=universe.read("fileIO.txt")
   print u.n, u.t, u.G, u.parts[0]



if __name__ == "__main__":
   main()



