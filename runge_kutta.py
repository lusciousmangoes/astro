import numpy as np


def RungeKutta(r_start, r_end, h, P0, M0, L0, T0, rho0, P_func, M_func, L_func, T_func, stop, data):
	
    r = np.arange(r_start,r_end,h)
    P = np.arange(r_start,r_end,h)
    M = np.arange(r_start,r_end,h)
    L = np.arange(r_start,r_end,h)
    T = np.arange(r_start,r_end,h)


    P[0] = P0
    M[0] = M0
    L[0] = L0
    T[0] = T0

		
    for i in range(0,len(r)-1):

        p1 = P_func(r[i], T[i], M[i], L[i])*h 
        m1 = M_func(r[i], T[i], P[i], L[i])*h 
        l1 = L_func(r[i], T[i], M[i], P[i])*h 
        t1 = T_func(r[i], P[i], M[i], L[i])*h 

        p2 = P_func(r[i]+h/2.0, T[i]+t1/2.0, M[i]+m1/2.0, L[i]+l1/2.0) * h 
        m2 = M_func(r[i]+h/2.0, T[i]+t1/2.0, P[i]+p1/2.0, L[i]+l1/2.0) * h 
        l2 = L_func(r[i]+h/2.0, T[i]+t1/2.0, M[i]+m1/2.0, P[i]+p1/2.0) * h 
        t2 = T_func(r[i]+h/2.0, P[i]+p1/2.0, M[i]+m1/2.0, L[i]+l1/2.0) * h 

        p3 = P_func(r[i]+h/2.0, T[i]+t2/2.0, M[i]+m2/2.0, L[i]+l2/2.0) * h 
        m3 = M_func(r[i]+h/2.0, T[i]+t2/2.0, P[i]+p2/2.0, L[i]+l2/2.0) * h 
        l3 = L_func(r[i]+h/2.0, T[i]+t2/2.0, M[i]+m2/2.0, P[i]+p2/2.0) * h 
        t3 = T_func(r[i]+h/2.0, P[i]+p2/2.0, M[i]+m2/2.0, L[i]+l2/2.0) * h

        p4 = P_func(r[i]+h, T[i]+t3, M[i]+m3, L[i]+l3) * h 
        m4 = M_func(r[i]+h, T[i]+t3, P[i]+p3, L[i]+l3) * h 
        l4 = L_func(r[i]+h, T[i]+t3, M[i]+m3, P[i]+p3) * h 
        t4 = T_func(r[i]+h, P[i]+p3, M[i]+m3, L[i]+l3) * h

        P[i+1] = P[i] + p1/6.0 + p2/3.0 + p3/3.0 + p4/6.0
        M[i+1] = M[i] + m1/6.0 + m2/3.0 + m3/3.0 + m4/6.0
        L[i+1] = L[i] + l1/6.0 + l2/3.0 + l3/3.0 + l4/6.0
        T[i+1] = T[i] + t1/6.0 + t2/3.0 + t3/3.0 + t4/6.0

        if stop(P[i+1],M[i+1],L[i+1],T[i+1],data):
            return P[0:i],M[0:i],L[0:i],T[0:i]

	#print("It should not reach this point")
    return P,M,L,T
