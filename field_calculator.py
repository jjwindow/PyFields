"""
PyFields field_calculator.py

13.10.2020

This file contains the explicit field strength calculator using the spherical harmonic expansion
model of planetary magnetic fields, in the form given in Connerney (1993).
"""

import numpy as np
import scipy as sp

###################### GLOBAL DEFINITIONS ############################

# Legendre (n,m) functions
lgd = np.array([[lambda theta: 1, lambda theta: 0, lambda theta: 0, lambda theta: 0], 
                    [lambda theta: np.cos(theta), lambda theta: np.sin(theta), lambda theta: 0, lambda theta: 0], 
                    [lambda theta: (3/2)*((np.cos(theta))**2-(1/3)), lambda theta: (3**0.5)*(np.cos(theta))*(np.sin(theta)), lambda theta: ((3**0.5)/2)*(np.sin(theta))**2, lambda theta: 0], 
                    [lambda theta: (5/2)*(np.cos(theta))*((np.cos(theta))**2 - (9/15)), lambda theta: ((5*(3**0.5))/(2**1.5))*(np.sin(theta))*((np.cos(theta))**2 - (3/15)), 
                    lambda theta: ((15**0.5)/2)*(np.cos(theta))*((np.sin(theta))**2), lambda theta: ((5**0.5)/(2**1.5))*((np.sin(theta))**3)]])

lgd_prime = np.array([[lambda theta: 0, lambda theta: 0, lambda theta:0],
                        [lambda theta: -np.sin(theta), lambda theta: np.cos(theta), lambda theta: 0],
                        [lambda theta: -(3/2)*np.sin(2*theta), lambda theta: (3**0.5)*((np.cos(theta))**2 - (np.sin(theta))**2), lambda theta: ((3**0.5)/2)*(np.sin(2*theta))]])

def legendre(n, m, th):
    """
    Function to return value of Legendre polynomial degree n, order m, calculated at angle theta.
    """
    return lgd[n][m](th)

def legendre_prime(n, m, th):
    """
    Function to return value of differentiated Legendre polynomial degree n, order m, calculated at angle theta.
    """
    return lgd_prime[n][m](th)

def B(p, field_coeffs):
    """
    Finds magnetic field strength at given (t, th, ph) co-ords for a given set of harmonic expansion 
    coefficients. Returns vector of components as a tuple.
    """

    # # Commented out with deprecation of planet kwarg.
    # # Set args tuple according to planet kwarg
    # if planet == "Uranus":
    #     args = (r, th, ph, a_U, g_U, h_U)
    # elif planet == "Neptune":
    #     args = (r, th, ph, a_N, g_N, h_N)
    # else:
    #     raise Exception("Keyword argument 'planet' must be 'Uranus' or 'Neptune'.")
    
    r, th, ph = p[0], p[1], p[2]
    args = (r, th, ph, *field_coeffs)
    # breakpoint()
    
    # Field component functions
    def B_rad(r, th, ph, a, g, h):
        """
        Radial magnetic field component. Formula from Connerney (1993).
        """
        B_rad_result= .0
        for n in range(0,3):
            for m in range(0, n+1):
                B_rad_result += (n+1)*((a/r)**(n+1))*(g[n][m]*np.cos(m*ph) + h[n][m]*np.sin(m*ph))*legendre(n, m, th)
        return B_rad_result

    def B_theta(r, th, ph, a, g, h):
        """
        Latitudinal magnetic field component. Formula from Connerney (1993).
        """
        B_theta_result= .0
        for n in range(0,3):
            for m in range(0, n+1):
                B_theta_result += -(a/r)**(n+2)*(g[n][m]*np.cos(m*ph) + h[n][m]*np.sin(m*ph))*legendre_prime(n, m, th)
        return B_theta_result

    def B_phi(r, th, ph, a, g, h):
        """
        Longitudinal magnetic field component. Formula from Connerney (1993).
        """
        B_phi_result= .0
        for n in range(0,3):
            for m in range(0, n+1):
                B_phi_result += (1/(np.sin(th)))*m*(a/r)**(n+2)*(g[n][m]*np.sin(m*ph) - h[n][m]*np.cos(m*ph))*legendre(n, m, th)
        return B_phi_result

    return np.array([B_rad(*args), B_theta(*args), B_phi(*args)])


# TESTING

