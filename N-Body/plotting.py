import numpy as np
import matplotlib.pyplot as plt

(t,E,me,ve,ea,ma,ju,sa,ur,ne)=np.loadtxt("data_orbits.dat",unpack=True)

#plt.plot(rs,psis,linewidth=3)
#plt.gcf().subplots_adjust(bottom=0.15)

#plt.xlabel(r'$x$',fontsize=24)
#plt.ylabel(r'$xe^{-x^2}$',fontsize=24)
#plt.xticks([0,2.1*10**-15],['0','a'])
#plt.yticks([])



plt.plot(t,me,linewidth=1,color = 'r')
plt.plot(t,ve,linewidth=1,color = 'b')
plt.plot(t,ea,linewidth=1,color = 'b')
plt.plot(t,ma,linewidth=1,color = 'b')
plt.plot(t,ju,linewidth=1,color = 'b')
plt.plot(t,sa,linewidth=1,color = 'b')
plt.plot(t,ur,linewidth=1,color = 'b')
plt.plot(t,ne,linewidth=1,color = 'b')
plt.savefig("orbits_plot.pdf")


