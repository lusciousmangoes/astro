#Defines all global variables
import scipy.constants

def init(x, y, z):
    global a, gamma, m_H, k, G, c, f_pp, f_3a, f_ff, mu, X, Y, Z
    a = 7.565767 * 10 **(-16)
    gamma = 5./3.
    m_H = scipy.constants.m_p
    k = scipy.constants.k
    G = scipy.constants.G
    c = scipy.constants.c
    f_pp = 1.
    f_3a = 1.
    g_ff = 1.
    mu = (2*x + 3./4. * y + 1./2. * z)**-1
    X = x
    Y = y
    Z = z    
