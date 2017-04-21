from numpy import *
import matplotlib.pyplot as plt


def read_groups_catalogue(filename):
  """
  Read the "fof_special_catalogue" files and return 4 arrays:
  
  GroupLen	: the size of each group
  GroupOffset	: the offset of the first particle in each group
                  as it will be found in the file fof_special_particles
  GroupMass	: the mass of each group (in code unit)
  GroupCM	: the center of mass of each group
  """


  f = open(filename,'r')

  Ngroups = fromstring(f.read(4),int32)[0]

  GroupLen    = fromstring(f.read(4*Ngroups),int32)
  GroupOffset = fromstring(f.read(4*Ngroups),int32)
  GroupMass   = fromstring(f.read(4*Ngroups),float32)

  GroupCM     = fromstring(f.read(3*4*Ngroups),float32)
  GroupCM.shape  = (Ngroups,3)


  GroupNspecies = fromstring(f.read(3*4*Ngroups),int32)
  GroupNspecies.shape  = (Ngroups,3)

  GroupMspecies = fromstring(f.read(3*4*Ngroups),float32)
  GroupMspecies.shape  = (Ngroups,3)

  GroupSfr   = fromstring(f.read(4*Ngroups),float32)

  GroupMetallicities = fromstring(f.read(2*4*Ngroups),float32)
  GroupMetallicities.shape  = (Ngroups,2)

  Mcold = fromstring(f.read(4*Ngroups),float32) 

  SigmaStars= fromstring(f.read(4*Ngroups),float32) 
  SigmaDM= fromstring(f.read(4*Ngroups),float32) 

  f.close()

  return GroupLen,GroupOffset,GroupMass,GroupCM
  
  
  
  
def read_groups_particles(filename):
  """
  Read the "fof_special_particles" files and return
  an array of the positions of each particles belonging
  to a group.
  """
  
  f = open(filename,'r')

  Ntot = fromstring(f.read(4),int32)[0]
  Pos	  = fromstring(f.read(3*4*Ntot),float32)
  Pos.shape  = (Ntot,3)
  f.close()
  
  return Pos





filename="fof_special_catalogue_036"

length,offset,mass,cm = read_groups_catalogue(filename)
#print mass,cm


plt.xlim(0,50)
plt.ylim(0,50)
x=[]
y=[]
for a in cm:
    x.append(a[0])
    y.append(a[1])
plt.scatter(x,y,s=length/2.)
#plt.show()
plt.savefig("N64.pdf")

