from numpy import exp, power, logspace, log10, zeros, sqrt, pi
import matplotlib.pyplot as plt
from matplotlib import rcParams

plt.rc('font', family='serif')
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
rcParams.update({'figure.autolayout': True})

#constants
G = 6.673e-11
c = 2.99792458e08
h = 6.62606876e-34
k_B = 1.3806503e-23
sigma = 2.0 * power(pi, 5.0) * power(k_B, 4.0) / (15.0 * c * c * power(h, 3.0))
a_rad = 4.0 * sigma / c
m_e = 9.10938188e-31
m_p = 1.67262158e-27
m_n = 1.67492716e-27
m_H = 1.673532499e-27
M_sun = 1.9891e30
L_sun = 3.828e26
R_sun = 6.95508e8
Te_sun = power(L_sun / (4.0 * pi * R_sun * R_sun * sigma), 0.25)
gamma = 5.0 / 3.0


# mean molecular weight, assuming complete ionization
def molecular_weight(X, Y, Z):
    return 1.0 / ( 2.0 * X + 0.75 * Y + 0.5 * Z)

# Nuclear reaction energy generation formulae
def epsilon(rho, T, X, Y, Z):
    T_6 = T / 1.0e6
    f_pp = 1.0
    psi_pp = 1.0 + 1.412e8 * (1.0/X - 1.0) * exp(-49.98 * power(T_6, -1.0/3.0))
    C_pp = 1.0 + 0.0123 * power(T_6, 1.0/3.0) + 0.0109 * power(T_6, 2.0/3.0) + 0.000938 * T_6
    
    e_pp = 0.241 * rho * X * X * f_pp * psi_pp * C_pp * power(T_6, -2.0/3.0) * exp(-33.80 * power(T_6, -1.0/3.0))
    
    C_cno = 1.0 + 0.0027 * power(T_6, 1.0/3.0) - 0.00778 * power(T_6, 2.0/3.0) - 0.000149 * T_6
    X_cno = 0.5 * Z
    
    e_cno = 8.67e20 * rho * X * X_cno * C_cno * power(T_6, -2.0/3.0) * exp(-152.28 * power(T_6, -1.0/3.0))
    
    T_8 = T/1.0e8
    f_3a = 1.0
    e_3a = 50.9 * rho * rho * Y * Y * Y / (T_8 * T_8 * T_8) * f_3a * exp(-44.027 / T_8)
    
    #return e_pp, e_cno, e_3a
    return e_pp + e_cno + e_3a
    

# opacity    
def opacity(rho, T, X, Y, Z):
    
    t_by_g_bf = 0.708 * power(rho * (1.0 + X), 0.2)
    kappa_bf = 4.34e21  / t_by_g_bf * Z * (1.0 + X) * rho * power(T, -3.5)
    
    g_ff = 1.0
    kappa_ff = 3.68e18 * g_ff * (1.0 - Z) * (1.0 + X) * rho * power(T, -3.5)
    
    kappa_es = 0.02 * (1.0 + X)
    
    kappa_Hm = 0.0
    if T > 3000.0 and T < 6000.0 and rho > 1e-7 and rho < 1e-2 and Z > 0.001 and Z < 0.003:
        kappa_Hm = 7.9e-34 * (Z / 0.02) * sqrt(rho) * power(T, 9.0)
    
    #return kappa_bf, kappa_ff, kappa_es, kappa_Hm
    return (kappa_bf + kappa_ff + kappa_es + kappa_Hm)

# density from equation of state    
def density(P, T, X, Y, Z):
    P_rad = a_rad * T * T * T * T / 3.0
    P_gas = P - P_rad
    mu = molecular_weight(X, Y, Z)
    rho = P_gas * mu * m_H / k_B / T

    return rho    




# test functions

def test_epsilon():
    rho = 162200
    N = 100
    T = logspace(6, 9, N)
    
    e_pp = zeros(N)
    e_cno = zeros(N)
    e_3a = zeros(N)
    e_tot = zeros(N)
    for i in range(N):
        e_pp[i], e_cno[i], e_3a[i] = epsilon(rho, T[i], 0.73, 0.25, 0.02)
        e_tot[i] = e_pp[i] + e_cno[i] + e_3a[i]
        if e_3a[i] == 0.0:
            e_3a[i] = 1e-300
        
    fig = plt.figure(figsize=(4, 3), dpi=200)
    ax = fig.add_subplot(1,1,1)
    ax.loglog(T, e_pp, color="black")
    ax.loglog(T, e_cno, color="red")
    ax.loglog(T, e_3a, color="blue")
    #ax.loglog(T, e_tot, color="black", lw="2")
    ax.set_ylim(1e-10, 1e14)
    ax.set_xlabel(r'$T$')
    ax.set_ylabel(r'$\epsilon$')    
    plt.savefig("epsilon.pdf")
    
def test_opacity():
    rho = 1408
    T = 1e4
    
    print opacity(rho, T, 0.73, 0.25, 0.02)
    

def main():
    T0 = 1.57e7
    P0 = 2.342e16
    
    for i in range(100):
        T = T0 - i/99.0 * T0
        P = P0 - i/99.0 * P0
        rho = density(P, T, 0.73, 0.25, 0.02)
        kappa = opacity(rho, T, 0.73, 0.25, 0.02)
        
        print T, P, rho, kappa

if __name__ == "__main__":
    main()
