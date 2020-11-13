"""
PyFields all_funcs.py

01.11.2020

All modules in PyFields combined into one file to improve execution time.
"""

import numpy as np
import numba
from tqdm import tqdm
import matplotlib.pyplot as plt

######################### GLOBAL DEFINITIONS #############################

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

########################## FIELD CALCULATOR ###############################
# Field component functions
@numba.njit
def _B_rad_(r, th, ph, a, g, h):
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
def _B_theta_(r, th, ph, a, g, h):
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
def _B_phi_(r, th, ph, a, g, h):
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
    out_array = np.array([_B_rad_(*args), _B_theta_(*args), _B_phi_(*args)])
    return out_array

############################## STEPPER #################################
@numba.njit
def B_mag(B, r, th):
    # r = p[0]
    # th = p[1]
    magnitude = np.sqrt(B[0]**2 + (r * B[1])**2 + (r * np.sin(th) * B[2])**2)
    # magnitude = B[0]
    return magnitude

def RK4(p_0, B_0, ds, field_coeffs):
    """
    Given starting coordinates r, th, ph, performs an RK4 step of size ds to
    get to follow the field to a new postion vector.

    ## ALL ISNTANCES OF PLANET -> field_coeffs ##
    """
    r, th, ph = p_0[0], p_0[1], p_0[2]
    # Field vector at starting point
    # take unit vector
    v_0 = B_0/B_mag(B_0, p_0[0], p_0[1])
    
    # First Euler step
    p_1 = p_0 + 0.5*ds*v_0
    B_1 = B(p_1, field_coeffs)
    v_1 = B_1/B_mag(B_1, p_1[0], p_1[1])

    # First correction step
    p_2 = p_0 + 0.5*ds*v_1
    B_2 = B(p_2, field_coeffs)
    v_2 = B_2/B_mag(B_2, p_2[0], p_2[1])

    # Second correction step
    p_3 = p_0 + ds*v_2
    B_3 = B(p_3, field_coeffs)
    v_3 = B_3/B_mag(B_3, p_3[0], p_3[1])

    p_next = p_0 + ds*(v_0 + 2*v_1 + 2*v_2 + v_3)/6
    B_next = B(p_next, field_coeffs)

    return p_next, B_next
    
def stepsizer():
    """
    We need a function that can calculate the field gradient at a point. We need to figure out
    what inputs are required for this. I think it should be possible to calculate the gradient
    from the vectors already given in RK4(). In which case, it should be done in RK4 and passed as an 
    argument here, where this function will return a stepsize inversely proportional to the grad mag.

    If it's more complicated, this can become the general Euler double-step which focuses on precision.
    This is more general but I think more computationally expensive.
    """
    pass

################################### WRAPPER ###################################

def field_trace(start_pos, field_coeffs, ds, max_iter, axes = "Cartesian"):
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
    p_arr = np.asarray([p for p in p_arr if np.any(p)])
    B_arr = np.asarray([b for b in B_arr if np.any(b)])
    if axes == "Cartesian":
        x, y = map(list, zip(*[(r*np.sin(theta), r*np.cos(theta)) for r, theta in zip(p_arr[:, 0], p_arr[:, 1])]))
        return x, y
    else:
        return p_arr, B_arr

def multilines(num, th_min = 0, th_max = np.pi, coeffs = dipole, ds = 0.01, maxits = 100000):
    """
    Plots 'num' (int) field lines for equally spaced theta values between th_min and th_max.
    Field lines calculated using field coefficients given by coeffs (tuple), stepsize ds (float),
    and terminating after maxits (int). Use plt.show() to display plot after calling.

    MOVE TO ALL_FUNCS EVENTUALLY
    """
    th_values = np.linspace(th_min, th_max, num)
    with tqdm(total = len(th_values), desc=f"THETA {round(th_min/np.pi, 2)}*pi TO {round(th_max/np.pi, 2)}*pi") as bar:
        for th in th_values:
            if th==0 or th==np.pi or th==2*np.pi:
                pass
            else:
                x, y = field_trace([1., th, 0.], coeffs, ds, maxits)
                if y[0] > y[-1]:
                    colour = 'r'
                else:
                    colour = 'b'
                plt.plot(x, y, color = colour)
                #  print(th)
            bar.update()