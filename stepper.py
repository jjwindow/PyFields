"""
PyFields stepper.py
20.10.20

Runge-Kutte IV field line tracing function. Uses B() function in field_calculator.py to
calculate the magnetic field vector at every step.
"""

import numpy as np
from field_calculator import *

def B_mag(B, p):
    return np.sqrt(B[0]**2 + p[0]**2 * B[1]**2 + (p[0]*np.sin(p[1]))**2 * B[2])

def RK4(p_0, ds, planet):
    """
    Given starting coordinates r, th, ph, performs an RK4 step of size ds to
    get to follow the field to a new postion vector.
    """
    # Field vector at starting point
    B_0 = B(r, th, ph, planet)
    # take unit vector
    v_0 = B_0/B_mag(B_0, p_0)
    
    # First Euler step
    p_1 = p_0 + 0.5*ds*v_0
    B_1 = B(*p_1, planet = planet)
    v_1 = B_1/B_mag(B_1, p_1)

    # First correction step
    p_2 = p_0 + 0.5*ds*v1
    B_2 = B(*p_2, planet=planet)
    v_2 = B_2/B_mag(B_2, p_2)

    # Second correction step
    p_3 = p_0 + ds*v_2
    B_3 = B(*p_3, planet=planet)
    v_3 = B_3/B_mag(*p_3, planet=planet)

    p_next = p_0 + ds*(v_0 + 2*v_1 + 2*v_2 + v_3)/6

    return p_next
    

