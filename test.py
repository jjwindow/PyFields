"""
PyFields test.py

13.10.2020

Repository placeholder for field tracing module.
"""

import numpy as np
import scipy as sp

#th=0
lgd=np.array([[1, 0, 0, 0], [np.cos(th), np.sin(th), 0, 0], [(3/2)*((np.cos(th))**2-(1/3)), (3**0.5)*(np.cos(th))*(np.sin(th)), ((3**0.5)/2)*(np.sin(th))**2, 0], [(5/2)*(np.cos(th))*((np.cos(th))**2 - (9/15)), ((5*(3**0.5))/(2**1.5))*(np.sin(th))*((np.cos(th))**2 - (3/15)), ((15**0.5)/2)*(np.cos(th))*((np.sin(th))**2), ((5**0.5)/(2**1.5))*((np.sin(th))**3)]], dtype='float64')
#print(lgd)

"""made these arrays using data from for Uranus"""
g=np.array([[0., 0., 0.], [0.11893, 0.11579, 0], [-0.06030, -0.12587, 0.00196]], dtype='float64')
h=np.array([[0., 0., 0.], [0., -0.15648, 0.], [0., 0.06116, 0.04759]], dtype='float64')

def B_rad(a, r, th, ph, g, h):
    B_rad_result= .0
    for n in range(0,3):
        for m in range(0, n+1):
            B_rad_result += (n+1)*((a/r)**(n+1))*(g[n][m]*np.cos(m*ph) + h[n][m]*np.sin(m*ph))*lgd[n][m]
    return B_rad_result

def B_theta(a, r, th, ph, g, h):
    B_theta_result= .0
    for n in range(0,3):
        for m in range(0, n+1):
            B_theta_result += -(a/r)**(n+2)*(g[n][m]*np.cos(m*ph) + h[n][m]*np.sin(m*ph))*lgd[n][m]
    return B_theta_result

def B_phi(a, r, th, ph, g, h):
    B_phi_result= .0
    for n in range(0,3):
        for m in range(0, n+1):
            B_phi_result += (1/(np.sin(th)))*m*(a/r)**(n+2)*(g[n][m]*np.sin(m*ph) - h[n][m]*np.cos(m*ph))*lgd[n][m]
    return B_phi_result

"""just to let you know, B_theta is wrong as i need to take the derivative of lgd but just trying to figure out how to do that, also i
think the whole thing may be wrong lol as i am getting an error"""