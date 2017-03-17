import numpy as np
import matplotlib.pyplot as plt


(r,m,t,l)=np.loadtxt("stars.dat",unpack=True)


l=l/3.848/10.**26.

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.invert_xaxis()
ax1.set_xlabel(r'$T\ (K)$',fontsize=17)
ax1.set_ylabel(r'$L\ (L_\odot)$',fontsize=17)

cm = plt.cm.coolwarm

h = plt.scatter(t, l, marker='o', facecolors=m/1.99/10.**30., edgecolors='k', s=r/7.**8., label='the data',cmap=cm)
plt.grid()
plt.gray()



cbar = plt.colorbar(mappable=h, ax=ax1)
cbar.set_label(r'$M\:(M_\odot)$',fontsize=17)

plt.savefig('extra_hr.pdf')
