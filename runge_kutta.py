def RungeKutta(r_start,r_end,h,P0,T0,P_func,M_func,L_func,T_func,stop,data):
	
	r = np.arange(r_start,r_end,h)
	P = np.arange(r_start,r_end,h)
	M = np.arange(r_start,r_end,h)
	L = np.arange(r_start,r_end,h)
	T = np.arange(r_start,r_end,h)

	P[0] = P0
	M[0] = 0
	L[0] = 0
	T[0] = T0
	
	for i in range(0,len(r)-1):

	    k1 = P(r[i],P[i],T[i])*h 
	    l1 = M(r[i],y[i],z[i])*h

	    k2 = P(r[i]+h/2.0,P[i]+k1/2.0,z[i]+l1/2.0)*h 
	    l2 = M(r[i]+h/2.0,y[i]+k1/2.0,z[i]+l1/2.0)*h

	    k3 = P(r[i]+h/2.0,y[i]+k2/2.0,z[i]+l2/2.0)*h 
	    l3 = M(r[i]+h/2.0,y[i]+k2/2.0,z[i]+l2/2.0)*h

	    k4 = P(r[i]+h,y[i]+k3,z[i]+l3)*h 
	    l4 = M(r[i]+h,y[i]+k3,z[i]+l3)*h

	    z[i+1] = z[i] + l1/6.0 + l2/3.0 + l3/3.0 + l4/6.0
	    y[i+1] = y[i] + k1/6.0 + k2/3.0 + k3/3.0 + k4/6.0

	    if stop(x[i+1],y[i+1],z[i+1],y[0]):
	    	return x[0:i],y[0:i]

	#print("It should not reach this point")
	return x,y
