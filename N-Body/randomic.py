from nbody import *
from random import *

n=10000.
part=[]
i=0

while i<n :
    a = uniform(-1,1)
    b = uniform(-1,1)
    c = uniform(-1,1)
    if a**2+b**2+c**2 <= 1.0:
        part.append(particle(1./n,[a,b,c],[0.,0.,0.]))
        i+=1

u=universe(n,0.,part)

u.writebinary('randomic.dat')

