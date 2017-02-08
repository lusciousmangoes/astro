import astro_4910 as phy

#euler(x_start,x_end,h,y0,z0,f,g,stop,data)

start = 0.0001
end = 10.
h = 0.01
y0 = 1.
z0 = 0.


def f(x,y,z):
   ret_val = z
   return ret_val


def g(x,y,z):
   n = 1
   ret_val = -y**n - 2.*z/x
   return ret_val

def stop(x,y,z,data):
   if y < h:
      return True
   else:
      return False



x,y,z=phy.runge_kutta(start,end,h,y0,z0,f,g,stop,1)

phy.plot(x,y,"x","y","")
#phy.plot(x,y,"x","y","example_fig_name.pdf")
