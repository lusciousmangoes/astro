import numpy as np
import matplotlib.pyplot as plt


(r,m,t,l)=np.loadtxt("stars_positive_luminosities.dat",unpack=True)


l=l/3.848/10.**26.

fig = plt.figure()
ax1 = fig.add_subplot(111)
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.invert_xaxis()
ax1.set_xlabel(r'$T\ (K)$',fontsize=17)
ax1.set_ylabel(r'$L\ (L_\odot)$',fontsize=17)



ax1.scatter(t, l, marker='o', facecolors='none', edgecolors='k', s=r/8.**8., label='the data')
plt.grid()
plt.savefig('new_bad_hr.pdf')
