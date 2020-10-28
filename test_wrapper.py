"""
PyFields test_wrapper.py

26.10.2020

First iteration of a wrapper function to test the stepper and field_calculator modules.
"""
import numpy as np
from stepper import *
from field_calculator import *

# Uranus Coefficients
g_U = np.array([[0., 0., 0.], [0.11893, 0.11579, 0.], [-0.06030, -0.12587, 0.00196]])
h_U = np.array([[0., 0., 0.], [0., -0.15648, 0.], [0., 0.06116, 0.04759]])
a_U = 1

uranus = (a_U, g_U, h_U)

# Neptune Coefficients
g_N = np.array([[0., 0., 0., 0.], [0.09732, 0.03220, 0., 0.], [0.07448, 0.00664, 0.04499, 0.], [-0.06592, 0.04098, -0.03581, 0.00484]])
h_N = np.array([[0., 0., 0., 0.], [0., 0.09889, 0., 0.], [0., 0.11230, -0.00070, 0.], [0., -0.03669, 0.01791, -0.00770]])
a_N = 1

neptune = (a_N, g_N, h_N)

# Dipole coefficients
g_D = np.array([[0., 0., 0., 0.], [1., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.]])
h_D = np.array([[0., 0., 0., 0.], [1., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.]])
a_D = 1

dipole = (a_D, g_D, h_D)

def field_trace(start_pos, field_coeffs, ds, max_iter):
    """
    Function to trace a field line given a starting positon.

    PARAMS
    ---------------------------------------------------------------------------------------------------
    start_pos       -   list or np.array; starting position in spherical coordinates from which to 
                        trace the field.
    field_coeffs    -   tuple of lists or np.arrays (2D);  2D array of coefficients for the spherical 
                        harmonic expansion. Must be of form (a, g, h).
    ds              -   float; stepsize to trace the field with, taken as a constant for all steps 
                        currently with intention of adding variable step-size to the tracer.
    max_iter        -   int; maximum number of iterations.

    RETURNS
    ----------------------------------------------------------------------------------------------------
    p_arr           -   np.array (2D); array of (r, th, ph) coordinates which fall on the traced field
                        line.
    B_arr           -   np.array (2D); array of field vectors for every point in p_arr.
    """

    B_0 = B(start_pos, field_coeffs)
    p_0 = start_pos

    p_arr, B_arr = np.asarray([np.zeros(3) for _ in range(max_iter)]), np.asarray([np.zeros(3) for _ in range(max_iter)])
    p_arr[0] = p_0
    B_arr[0] = B_0

    it = 1
    while (p_0[0] >= 1.) and (it < max_iter):
        p_next, B_next = RK4(p_0, B_0, ds, field_coeffs)
        p_arr[it] = p_next
        B_arr[it] = B_next
        p_0, B_0 = p_next, B_next
        it += 1

    #return B_arr
    x, y = map(list, zip(*[(r*np.sin(theta), r*np.cos(theta)) for r, theta in zip(p_arr[:, 0], p_arr[:, 1])]))
    return x, y