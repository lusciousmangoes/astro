def RungeKutta(r_start,r_end,h,P0,T0,f,g,h,l,stop,data):
	
	r = np.arange(x_start,x_end,h)
	P = np.arange(x_start,x_end,h)
	M = np.arange(x_start,x_end,h)
	L = np.arange(x_start,x_end,h)
	T = np.arange(x_start,x_end,h)

	P[0] = P0
	M[0] = 0
	L[0] = 0
	T[0] = T0
	
	for i in range(0,len(x)-1):

	    k1 = f(x[i],y[i],z[i])*h 
	    l1 = g(x[i],y[i],z[i])*h

	    k2 = f(x[i]+h/2.0,y[i]+k1/2.0,z[i]+l1/2.0)*h 
	    l2 = g(x[i]+h/2.0,y[i]+k1/2.0,z[i]+l1/2.0)*h

	    k3 = f(x[i]+h/2.0,y[i]+k2/2.0,z[i]+l2/2.0)*h 
	    l3 = g(x[i]+h/2.0,y[i]+k2/2.0,z[i]+l2/2.0)*h

	    k4 = f(x[i]+h,y[i]+k3,z[i]+l3)*h 
	    l4 = g(x[i]+h,y[i]+k3,z[i]+l3)*h

	    z[i+1] = z[i] + l1/6.0 + l2/3.0 + l3/3.0 + l4/6.0
	    y[i+1] = y[i] + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0

	    if stop(x[i+1],y[i+1],z[i+1],y[0]):
	    	return x[0:i],y[0:i]

	#print("It should not reach this point")
	return x,y
