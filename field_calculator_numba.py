"""
PyFields field_calculator.py

13.10.2020

This file contains the explicit field strength calculator using the spherical harmonic expansion
model of planetary magnetic fields, in the form given in Connerney (1993).

# MODIFIED FOR NUMBA APPLICATION #
"""

import numpy as np
import numba

###################### GLOBAL DEFINITIONS ############################

# Legendre (n,m) functions

# def legendre(n, m, th):
#     """
#     Function to return value of Legendre polynomial degree n, order m, calculated at angle theta.
#     """
#     f = lgd[n][m]
#     out = f(th)
#     return out

# def legendre_prime(n, m, th):
#     """
#     Function to return value of differentiated Legendre polynomial degree n, order m, calculated at angle theta.
#     """
#     f = lgd[n][m]
#     out = f(th)
#     return out

# Field component functions
@numba.njit
def B_rad(r, th, ph, a, g, h):
    """
    Radial magnetic field component. Formula from Connerney (1993).
    """
    lgd = [[1.,0.,0.,0.], 
            [np.cos(th),np.sin(th),0.,0.], 
            [(3/2)*((np.cos(th))**2-(1/3)),(3**0.5)*(np.cos(th))*(np.sin(th)),((3**0.5)/2)*(np.sin(th))**2,0.], 
            [(5/2)*(np.cos(th))*((np.cos(th))**2 - (9/15)),((5*(3**0.5))/(2**1.5))*(np.sin(th))*((np.cos(th))**2 - (3/15)), 
            ((15**0.5)/2)*(np.cos(th))*((np.sin(th))**2),((5**0.5)/(2**1.5))*((np.sin(th))**3)]]
    B_rad_result= .0
    for n in range(0,3):
        for m in range(0, n+1):
            B_rad_result += (n+1)*((a/r)**(n+1))*(g[n][m]*np.cos(m*ph) + h[n][m]*np.sin(m*ph))*lgd[n][m]
    return B_rad_result

@numba.njit
def B_theta(r, th, ph, a, g, h):
    """
    Latitudinal magnetic field component. Formula from Connerney (1993).
    """
    lgd_prime = [[0., 0., 0.],
                [-np.sin(th), np.cos(th), 0.],
                [-(3/2)*np.sin(2*th),(3**0.5)*((np.cos(th))**2 - (np.sin(th))**2),((3**0.5)/2)*(np.sin(2*th))]]

    B_theta_result= .0
    for n in range(0,3):
        for m in range(0, n+1):
            B_theta_result += -(a/r)**(n+2)*(g[n][m]*np.cos(m*ph) + h[n][m]*np.sin(m*ph))*lgd_prime[n][m]
    return B_theta_result

@numba.njit
def B_phi(r, th, ph, a, g, h):
    """
    Longitudinal magnetic field component. Formula from Connerney (1993).
    """
    lgd = [[1.,0.,0.,0.], 
            [np.cos(th),np.sin(th),0.,0.], 
            [(3/2)*((np.cos(th))**2-(1/3)),(3**0.5)*(np.cos(th))*(np.sin(th)),((3**0.5)/2)*(np.sin(th))**2,0.], 
            [(5/2)*(np.cos(th))*((np.cos(th))**2 - (9/15)),((5*(3**0.5))/(2**1.5))*(np.sin(th))*((np.cos(th))**2 - (3/15)), 
            ((15**0.5)/2)*(np.cos(th))*((np.sin(th))**2),((5**0.5)/(2**1.5))*((np.sin(th))**3)]]
    B_phi_result= .0
    for n in range(0,3):
        for m in range(0, n+1):
            B_phi_result += (1/(np.sin(th)))*m*(a/r)**(n+2)*(g[n][m]*np.sin(m*ph) - h[n][m]*np.cos(m*ph))*lgd[n][m]
    return B_phi_result

def B(p, field_coeffs):
    """
    Finds magnetic field strength at given (t, th, ph) co-ords for a given set of harmonic expansion 
    coefficients. Returns vector of components as a tuple.
    """
    r, th, ph = p[0], p[1], p[2]
    args = (r, th, ph, *field_coeffs)
    out_array = np.array([B_rad(*args), B_theta(*args), B_phi(*args)])
    return out_array


