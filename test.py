"""
PyFields test.py

13.10.2020

Repository placeholder for field tracing module.
"""

import numpy as np
import scipy as sp

#th=0
#lgd=np.array([[1, 0, 0, 0], [np.cos(th), np.sin(th), 0, 0], [(3/2)*((np.cos(th))**2-(1/3)), (3**0.5)*(np.cos(th))*(np.sin(th)), ((3**0.5)/2)*(np.sin(th))**2, 0], [(5/2)*(np.cos(th))*((np.cos(th))**2 - (9/15)), ((5*(3**0.5))/(2**1.5))*(np.sin(th))*((np.cos(th))**2 - (3/15)), ((15**0.5)/2)*(np.cos(th))*((np.sin(th))**2), ((5**0.5)/(2**1.5))*((np.sin(th))**3)]], dtype='float64')
#print(lgd)

###################### GLOBAL DEFINITIONS ############################

# Uranus Coefficients
g_U=np.array([[0., 0., 0.], [0.11893, 0.11579, 0], [-0.06030, -0.12587, 0.00196]], dtype='float64')
h_U=np.array([[0., 0., 0.], [0., -0.15648, 0.], [0., 0.06116, 0.04759]], dtype='float64')

# Legendre (n,m) functions
lgd = np.array([[lambda theta: 1, lambda theta: 0, lambda theta: 0, lambda theta: 0], 
                    [lambda theta: np.cos(theta), lambda theta: np.sin(theta), lambda theta: 0, lambda theta: 0], 
                    [lambda theta: (3/2)*((np.cos(theta))**2-(1/3)), lambda theta: (3**0.5)*(np.cos(theta))*(np.sin(theta)), lambda theta: ((3**0.5)/2)*(np.sin(theta))**2, lambda theta: 0], 
                    [lambda theta: (5/2)*(np.cos(theta))*((np.cos(theta))**2 - (9/15)), lambda theta: ((5*(3**0.5))/(2**1.5))*(np.sin(theta))*((np.cos(theta))**2 - (3/15)), 
                    lambda theta: ((15**0.5)/2)*(np.cos(theta))*((np.sin(theta))**2), lambda theta: ((5**0.5)/(2**1.5))*((np.sin(theta))**3)]])

def legendre(n, m, th):
    """
    Function to return value of Legendre polynomial degree n, order m, calculated at angle theta.
    """
    return lgd[n][m](th)

def B_rad(a, r, th, ph, g, h):
    """
    Radial magnetic field component. Formula from Connerney (1993).
    """
    B_rad_result= .0
    for n in range(0,3):
        for m in range(0, n+1):
            B_rad_result += (n+1)*((a/r)**(n+1))*(g[n][m]*np.cos(m*ph) + h[n][m]*np.sin(m*ph))*legendre(n, m, th)
    return B_rad_result

def B_theta(a, r, th, ph, g, h):
    """
    Latitudinal magnetic field component. Formula from Connerney (1993).
    """
    B_theta_result= .0
    for n in range(0,3):
        for m in range(0, n+1):
            B_theta_result += -(a/r)**(n+2)*(g[n][m]*np.cos(m*ph) + h[n][m]*np.sin(m*ph))*legendre(n, m, th)
    return B_theta_result

def B_phi(a, r, th, ph, g, h):
    """
    Longitudinal magnetic field component. Formula from Connerney (1993).
    """
    B_phi_result= .0
    for n in range(0,3):
        for m in range(0, n+1):
            B_phi_result += (1/(np.sin(th)))*m*(a/r)**(n+2)*(g[n][m]*np.sin(m*ph) - h[n][m]*np.cos(m*ph))*legendre(n, m, th)
    return B_phi_result

"""just to let you know, B_theta is wrong as i need to take the derivative of lgd but just trying to figure out how to do that, also i
think the whole thing may be wrong lol as i am getting an error"""

# TESTING

print(B_theta(1, 1.3, np.pi, 0, g_U, h_U))