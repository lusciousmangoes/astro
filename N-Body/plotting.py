import numpy as np
import matplotlib.pyplot as plt

(t,ex,ey,ez,sx,sy,sz,E)=np.loadtxt("data.dat",unpack=True)

#plt.plot(rs,psis,linewidth=3)
#plt.gcf().subplots_adjust(bottom=0.15)

#plt.xlabel(r'$x$',fontsize=24)
#plt.ylabel(r'$xe^{-x^2}$',fontsize=24)
#plt.xticks([0,2.1*10**-15],['0','a'])
#plt.yticks([])



plt.plot(ex,ey,linewidth=1,color = 'r')
plt.plot(sx,sy,linewidth=1,color = 'b')
plt.savefig("data_plot.pdf")

plt.close()
plt.plot(t,E,linewidth=1,color = 'r')
plt.savefig("data_plot2.pdf")