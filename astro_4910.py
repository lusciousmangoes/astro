import numpy as np
import matplotlib.pyplot as plt


def integrate(x,f):

   a = 0.0
   i = 0
   
   while i < len(x)-1:
      a += ((f[i]+f[i+1])/2.0)*(x[i+1]-x[i])
      i += 1

   return a



def plot(x,y,x_label,y_label,name):

   plt.plot(x,y,linewidth = 2, color = "blue",)
   plt.xlabel(x_label)
   plt.ylabel(y_label)
   plt.axhline(y=0, color='black')
   plt.axvline(x=0, color='black')
   if name == "":
      plt.show()
   else:
      plt.gcf().subplots_adjust(left=0.2)
      plt.savefig(name)
   plt.close()




def euler(x_start,x_end,h,y0,z0,f,g,stop,data):

   x = np.arange(x_start,x_end,h)
   y = np.arange(x_start,x_end,h)
   z = np.arange(x_start,x_end,h)

   y[0] = y0
   z[0] = z0


   for n in range(0,len(x)-1):

      z[n+1] = g(x[n],y[n],z[n])*h + z[n]
      y[n+1] = f(x[n],y[n],z[n])*h + y[n]

      if stop(x[n+1],y[n+1],z[n+1],data):
         return x[0:n],y[0:n],z[0:n]

   return x,y,z




def runge_kutta(x_start,x_end,h,y0,z0,f,g,stop,data):


   x = np.arange(x_start,x_end,h)
   y = np.arange(x_start,x_end,h)
   z = np.arange(x_start,x_end,h)

   y[0] = y0
   z[0] = z0



   for n in range(0,len(x)-1):



      k1 = f(x[n],y[n],z[n])*h
      l1 = g(x[n],y[n],z[n])*h

      k2 = f(x[n]+h/2.,y[n]+k1/2.,z[n]+l1/2.)*h
      l2 = g(x[n]+h/2.,y[n]+k1/2.,z[n]+l1/2.)*h

      k3 = f(x[n]+h/2.,y[n]+k2/2.,z[n]+l2/2.)*h
      l3 = g(x[n]+h/2.,y[n]+k2/2.,z[n]+l2/2.)*h

      k4 = f(x[n]+h,y[n]+k3,z[n]+l3/2.)*h
      l4 = g(x[n]+h,y[n]+k3,z[n]+l3/2.)*h

      y[n+1] = (k1+2*k2+2*k3+k4)/6. + y[n]

      z[n+1] = (l1+2*l2+2*l3+l4)/6. + z[n]


      if stop(x[n+1],y[n+1],z[n+1],y0):
         return x[0:n],y[0:n],z[0:n]

   return x,y,z


