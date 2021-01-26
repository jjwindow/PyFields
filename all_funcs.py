"""
PyFields all_funcs.py

01.11.2020

All modules in PyFields combined into one file to improve execution time.
"""

import numpy as np
import numba
from tqdm import tqdm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd
import warnings
import os.path
from palettable.wesanderson import Aquatic2_5, Cavalcanti_5
from numpy.linalg import norm
import matplotlib as mat

######################### GLOBAL DEFINITIONS #############################

# Uranus Coefficients
g_U = np.array([[0., 0., 0.], [0.11893, 0.11579, 0.], [-0.06030, -0.12587, 0.00196]])
h_U = np.array([[0., 0., 0.], [0., -0.15648, 0.], [0., 0.06116, 0.04759]])
g_U_err = np.array([[0., 0., 0.], [0.001, 0.003, 0.], [0.00550, 0.00610, 0.005]])
h_U_err = np.array([[0., 0., 0.], [0., 0.0017, 0.], [0., 0.00360, 0.00810]])
a_U = 1

uranus = (a_U, g_U, h_U)
uranus_uncert = (a_U, g_U_err, h_U_err)

# Neptune Coefficients
g_N = np.array([[0., 0., 0.], [0.09732, 0.03220, 0.], [0.07448, 0.00664, 0.04499]])
h_N = np.array([[0., 0., 0.], [0., -0.09889, 0.], [0., 0.11230, -0.00070]])
g_N_err = np.array([[0., 0., 0.], [0.002, 0.0036, 0.], [0.0113, 0.0112, 0.0084]])
h_N_err = np.array([[0., 0., 0.], [0., 0.0011, 0.], [0., 0.003, -0.0034]])
a_N = 1

neptune = (a_N, g_N, h_N)
neptune_uncert = (a_N, g_N_err, h_N_err)

# Dipole coefficients
g_D = np.array([[0., 0., 0., 0.], [1., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.]])
h_D = np.array([[0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.]])
a_D = 1

dipole = (a_D, g_D, h_D)

# Quadrupole coefficients
g_Q = np.array([[0., 0., 0., 0.], [0., 0., 0., 0.], [1., 0., 0., 0.], [0., 0., 0., 0.]])
h_Q = np.array([[0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.], [0., 0., 0., 0.]])
a_Q = 1

quadrupole = (a_Q, g_Q, h_Q)

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

def RK4(p_0, B_0, ds, field_coeffs, back = False):
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

    if not back:
        p_next = p_0 + ds*(v_0 + 2*v_1 + 2*v_2 + v_3)/6
    else:
        p_next = p_0 - ds*(v_0 + 2*v_1 + 2*v_2 + v_3)/6

    B_next = B(p_next, field_coeffs)

    return p_next, B_next


################################### WRAPPER ###################################

def field_trace(start_pos, field_coeffs, ds, max_iter, axes = "Cartesian", back = False):
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
        p_next, B_next = RK4(p_0, B_0, ds, field_coeffs, back)
        p_arr[it] = p_next
        B_arr[it] = B_next
        p_0, B_0 = p_next, B_next
        it += 1
        iter_flag = (it == max_iter)

    if (iter_flag):
        return None

    p_arr = np.asarray([p for p in p_arr if np.any(p)])[:-1]
    B_arr = np.asarray([b for b in B_arr if np.any(b)])[:-1]

    if ((len(p_arr) < 3) or iter_flag):
        return None
    else:
        if axes == "Cartesian":
            x, z, y = spherical2cartesian(p_arr)
            return x, y, z
        else:
            return p_arr, B_arr

def spherical2cartesian(p_arr):
    """
    Converts (2d) list of spherical coordinates to 3 (1d) lists of cartesian coordinates for
    use in plotting functions.
    PARAMS
    ---------------------------------------------------------------------------------------
    p_arr   -   list or numpy array (2D); list of [r, theta, phi] coordinates.
    RETURNS
    ---------------------------------------------------------------------------------------
    x, y, z -   lists (1D); list of x-, y-, and z-coordinates (respectively) in Cartesian axes
                corresponding to the same points in space defined by p_arr.
    """
    x, z, y = map(list, zip(*[(r*np.sin(theta)*np.cos(phi), r*np.cos(theta), r*np.sin(theta)*np.sin(phi)) for r, theta, phi in zip(p_arr[:, 0], p_arr[:, 1], p_arr[:, 2])]))
    return x, y, z

def multilines(phi, num, th_min = 0, th_max = 2*np.pi, coeffs = dipole, ds = 0.01, maxits = 100000, plot = True):
    """
    Plots 'num' (int) field lines for equally spaced theta values between th_min and th_max.
    Field lines calculated using field coefficients given by coeffs (tuple), stepsize ds (float),
    and terminating after maxits (int). Use plt.show() to display plot after calling.
    """
    th_values = np.linspace(th_min, th_max, num, endpoint=False)
    field_lines = []
    with tqdm(total = len(th_values), desc=f"THETA {round(th_min/np.pi, 2)}*pi TO {round(th_max/np.pi, 2)}*pi") as bar:
        for th in th_values:
            if th==0 or th==np.pi or th==2*np.pi:
                pass
            else:
                field_line = field_trace([1., th, float(phi)], coeffs, ds, maxits)
                if field_line is not None:
                    (x, y, z) = field_line
                    if plot:
                        if y[0] > y[-1]:
                            colour = 'r'
                        else:
                            colour = 'b'
                        plt.plot(x, y, color = colour)
                    else:
                        field_lines.append(field_line)
            bar.update()
    # field_lines = np.asarray(field_lines)
    return field_lines

def multiline_3D(num_th, phi_array, coeffs = dipole, ds = 0.01, maxits = 100000):
    fig=plt.figure()
    ax = plt.axes(projection = '3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    for phi in phi_array:
        field_lines = multilines(phi, num_th, th_min=-np.pi/2, th_max=np.pi/2, coeffs=coeffs, ds=ds, maxits=maxits, plot=False)
        for field_line in field_lines:
            (x, y, z) = field_line
            ax.plot3D(x, y, z, color='b')


##################### ANALYTIC COMPARISONS #######################

def _analytic_field_point_(th_i, th, field = 'dipole'):
    """
    Calculates (x,y) coordinate at th for a field line whith starting coordinate (t, th, ph) = (1, th, 0).
    Also returns rflag, which is True if r <= 1 and false otherwise. This is to terminate calculation.
    """
    if field == 'dipole':
        def x(th_i, th):
            return ((np.sin(th)**3)/np.sin(th_i)**2)
        def y(th_i, th):
            return ((np.sin(th)**2 * np.cos(th))/np.sin(th_i)**2)
    else:
        def x(th_i, th):
            return (np.sin(th_i)**2 * np.cos(th_i))**(-0.5) * np.sqrt(np.sin(th)**2 * np.cos(th)) * np.sin(th)
        def y(th_i, th):
            return (np.sin(th_i)**2 * np.cos(th_i))**(-0.5) * np.sqrt(np.sin(th)**2 * np.cos(th)) * np.cos(th)
    
    x, y = x(th_i, th), y(th_i, th)
    rflag = (round((x**2 + y**2), 6) < 1)        # Boolean flag - is radial coord < 1?
    # print(rflag)
    return x, y, rflag
    
def analytic_field_line(th_i, ds, field = 'dipole'):
    th_range = np.arange(th_i, 2*np.pi, step=ds)
    # th_i_range = np.array([th_i for _ in th_range])
    # x_y_coords = [(x, y) for x, y, rflag in [_analytic_field_point_(th_i, th) for th_i, th in zip(th_i_range, th_range)] if not rflag]
    x_y_coords = []
    j = 0
    rflag = False
    # breakpoint()
    while (not rflag) and (j < len(th_range)):
        x, y, rflag = _analytic_field_point_(th_i, th_range[j], field)
        x_y_coords.append((x, y))
        j += 1
    # breakpoint()
    return x_y_coords

def _analytic_field_plot(th_min, th_max, numlines, ds, field = 'dipole'):
    th_start = np.linspace(th_min, th_max, numlines, endpoint=False)
    for th_i in th_start:
        coords = analytic_field_line(th_i, ds, field)
        x_arr, y_arr = map(list, zip(*coords))
        plt.plot(x_arr, y_arr, '--', color = 'k')

################# COORDINATE TRANSFORM ########################

def cartesian2latlong(x, y, z):
    """
    Convert 3D Cartesian coordinates to latitude-longitudes for 
    2D projection plots.
    PARAMS
    -----------------------------------------------------------------------
    x, y, z -   float; coordinates in planet-centred Cartesian
                system. Axis of planetary rotation aligned along z-axis.
    RETURNS
    -----------------------------------------------------------------------
    lat, long - float;
    """
    # Convert lists to arrays for vectorisation.
    # Ignores floats and arrays.
    args = [x, y, z]
    for i, elem in enumerate(args):
        if isinstance(elem, list):
            args[i] = np.asarray(elem)
    
    [x, y, z] = args
    r = np.sqrt(x**2 + y**2 + z**2)
    lat = np.arcsin(z/r)*(180/(np.pi))
    longt = np.arctan2(y, x)*(180/(np.pi))

    return lat, longt

##################### MOON SELECTOR ###############################
df = pd.read_csv('satellite_properties.csv')
df.set_index('Name', inplace=True)

def moon_selector(moon, *args):
    """
    Returns desired parameters for a given moon.

    PARAMS
    ------------------------------------------------------------------
    Possible values of 'moon':'
    'Miranda' - 'Ariel' - 'Umbriel' - 'Titania' - 'Oberon' - 'Triton'
    (not case sensitive.)
    *args:
    Passing no args returns full properties dictionary for that moon.
    Otherwise, *args are keys to return properties for the moon. Each
    argument should be a string.
    Valid args are:
    'Parent'    -   'inc'      -    'R'    -    'a'     -   'T'
    (parent         (inclination,  (Radius,   (scaled     (orbital
    planet, str)     radians)           km)    radius)   time period)

                        'coeffs'    -    'uncert'
    Spherical harmonic coefficients for parent planet / associated 
    uncertainties.
    (tuple of form (a, g, h).)
                        'parent_day'
    length of a day on the parent planet in units of Earth days

    RETURNS
    -------------------------------------------------------------------
    out_dict    -   dict; of type {'arg' : arg_value, ...} for all 'arg' 
                    passed as arguments.

    Note - Invalid arguments do not raise an error but deploy a warning.
    """
    if not isinstance(moon, str):
        raise TypeError("Positional argument 'moon' must be of type string.")
    
    # df = pd.read_csv('satellite_properties.csv')
    # Select coefficients to use 
    # df.set_index('Name', inplace=True)
    moon = moon.lower()
    if moon not in df.index.values:
        raise ValueError("'moon' must be one of the 5 major Uranian moons or 'triton'.")
    
    moon_dict = df.loc[moon].to_dict()
    parent = moon_dict['Parent']
    coeffs = (parent == 'Uranus')*uranus + (parent == 'Neptune')*neptune
    coeffs_uncert = (parent == 'Uranus')*uranus_uncert + (parent == 'Neptune')*neptune_uncert
    moon_dict['coeffs'] = coeffs
    moon_dict['uncert'] = coeffs_uncert

    if not len(args):
        return moon_dict
    
    out_tup = tuple([moon_dict[arg] for arg in args if arg in moon_dict.keys()])

    bad_args = [arg for arg in args if not (arg in moon_dict.keys())]

    if len(bad_args):
        warnings.warn(f'The following arguments are not in satellite_proprties and were not returned:\n {bad_args}')

    return out_tup

######################### TIMER #################################

import time
def functimer(func, args, n):
    """
    Times a function n times, displays and returns the average time taken.
    'args' are arguments to pass to the function being timed and should be 
    a tuple to be unpacked.
    RETURNS
    ----------------------------------------------------------------------

    """
    t = 0
    for _ in range(n):
        t_0 = time.time()
        func(*args)
        t += time.time() - t_0
    mean = t/n
    print(f"{func.__name__} Time ({n} run avg):\n{mean}")
    return mean

######################## COEFFICIENT UNCERTAINTIES ######################'''''

def random_footpoints(n, moon, pos, trueTrace = False):
    """
    A function that generates random magnetic field footpoints within the bounds 
    of the uncertainties of the magnetic field coefficients. New spherical harmonic 
    expansion coefficients are calculated using a pseudorandom number generator, 
    and n fieldlines are traced using these for a given moon-planet system at a 
    fixed colatitude. The array of footpoints for all fieldlines is returned.

    PARAMS
    --------------------------------------------------------------------------------
    n       -   int; number of random fieldlines to calculate.
    moon    -   str; name of the moon to calculate footpoint uncertainties for.
    phi     -   float; value of colatitude at which to start the fieldline.
    trueTrace   -   bool; trace the fieldline using the accepted g, h coeffs.

    RETURNS
    ---------------------------------------------------------------------------------
    footpoints - list; list of tuples, where each tuple is (x, y, z) position of a 
                 footpoint of a fieldline calculated from the random coefficients.
    """

    (R, coeffs, uncert) = moon_selector(moon, 'a', 'coeffs', 'uncert')
    # start_pos = [R, np.pi/2, phi]
    (a, g, h) = coeffs
    (a, g_err, h_err) = uncert
    
    # Trace the accepted fieldline if desired
    if trueTrace:
        # x, y, z = field_trace(start_pos, (a, g, h), 0.005, 200000)
        x, y, z = field_trace(pos, coeffs, 0.005, 200000)
        trueFoot_f = (x[-1], y[-1], z[-1])
        x, y, z = field_trace(pos, coeffs, 0.005, 200000, back=True)
        trueFoot_b = (x[-1], y[-1], z[-1])
    # initialise footpoints array
    footpoints_f = [0. for _ in range(n)]
    footpoints_b = [0. for _ in range(n)]
    with tqdm(total=n, desc=":)") as bar:
        for k in range(n):
            g_new = np.zeros((3,3))
            h_new = np.zeros((3,3))

            for i in range(3):
                for j in range(3):
                    # Ignore null coefficients
                    if g[i][j] == 0.:
                        pass
                    else:
                        # Generate random num between -1 and 1
                        r_1 = (np.random.random()-0.5)*2
                        # Use random num as multiplier on uncertainty, add
                        # to coefficients
                        g_new[i][j] = g[i][j] + g_err[i][j]*r_1
                        # Repeat with different randnum for h coeffs
                        r_2 = (np.random.random() - 0.5)*2
                        h_new[i][j] = h[i][j] + h_err[i][j]*r_2
            
            coeffs = (a, g_new, h_new)
            # Trace fieldline with new set of coefficients
            x, y, z = field_trace(pos, coeffs, 0.005, 200000)
            # Take fieldline footpoint
            footpoints_f[k] = (x[-1], y[-1], z[-1])
            x, y, z = field_trace(pos, coeffs, 0.005, 200000, back=True)
            footpoints_b[k] = (x[-1], y[-1], z[-1])
            bar.update()

    if trueTrace:
        return footpoints_f, footpoints_b, trueFoot_f, trueFoot_b
    else:
        return footpoints_f, footpoints_b

######################## ORBIT CALCULATION ########################

def orbit(moon, num, num_orbits, relative = False):      #num_orbits is how many sidereal orbits #num gives num of points in one sidereal orbit
    """
    Function to generate coordinates of an orbital path of a given satellite around its parent.
    Can calculate orbits in the sidereal rest frame or in the planet's rest frame.

    PARAMS
    -----------------------------------------------------------------------------------
    moon        -   str; name of one of the 5 Uranian moons, or Triton.
    num         -   int; number of time segments to plot per orbit, i.e - time resolution.
    num_orbits  -   float or int; number of orbits to trace. Only makes a difference for 
                    inclined orbits with relative = True.
    relative    -   bool; if false, orbit calculated is in sidereal rest frame, i.e - no
                    consideration of planetary rotation. If true, then planetary rotation
                    is calculated and orbit given is the path seen from a frame co-rotating 
                    with the parent planet.
    RETURNS
    ------------------------------------------------------------------------------------
    orbital_points  -   numpy array; array containing num + 1 points in spherical 
                        coordinates, determining the orbital path. Each point is a list
                        length 3, [r, theta, phi].
    """
    # Collect moon parameters
    (R, coeffs, period_moon, period_plan, incl) = moon_selector(moon, 'a', 'coeffs', 'T', 'parent_day', 'inc')
    incl = (np.pi/180) * incl       # convert inclination to radians
    omega_moon = (2*np.pi)/period_moon      # period -> frequency
    omega_plan = (2*np.pi)/period_plan
    t_step = period_moon/num 
    n = int(num*num_orbits)     # number of points to plot - int() covers non-whole num_orbits.

    orbital_points= [0 for i in range(n+1)]     # initialise output list
    T_arr = [i*t_step for i in range(n+1)]

    for t in T_arr:
        # angular argument of satellite in the plane of its orbit, more correctly called the 'argument of latitude'.
        phi_moon_orbit = omega_moon * t     
        # from Adam's eqns:
        theta = np.arccos(np.cos(phi_moon_orbit)*np.sin(np.pi-incl))    
        phi_moon_eq = np.arctan2(-1*np.sin(phi_moon_orbit), np.cos(phi_moon_orbit)*np.cos(np.pi - incl))
        # phi_moon_eq is latitude coordinate in equatorial plane.
        if phi_moon_eq < 0:
            # handles negative arctan2 output
            phi_moon_eq += 2*np.pi
        if relative:
            # changes to planet rest frame
            phi = phi_moon_eq - omega_plan * t
        else:
            phi = phi_moon_eq 
        # append point to list
        pos = [R, theta, phi]
        orbital_points[i] = pos
    return np.array(orbital_points), np.array(T_arr)

###### Calculating mean angular error ######

def angular_deviation(footpoints_f_arr, trueFoot_f_arr, footpoints_b_arr, trueFoot_b_arr):

    mean_ang_dev_f = []
    mean_lat_dev_f = []
    mean_long_dev_f = []
    mean_ang_dev_b = []
    mean_lat_dev_b = []
    mean_long_dev_b = []

    for i, (pos, trueFoot) in enumerate(trueFoot_f_arr):
        trueLat, trueLongt = cartesian2latlong(*trueFoot)
        ang_dev_f = []
        lat_dev_f = []
        long_dev_f = []
        (pos, fp_arr) = footpoints_f_arr[i]
        for fp in fp_arr:
            # breakpoint()
            lat, longt = cartesian2latlong(*fp)
            ang_dev_f.append(np.arccos(np.dot(fp, trueFoot)/(norm(fp)*norm(trueFoot))))
            lat_dev_f.append(trueLat - lat)
            long_dev_f.append(trueLongt - longt)
        mean_ang_dev_f.append((pos, np.mean(ang_dev_f)))
        mean_lat_dev_f.append((pos, np.mean(lat_dev_f)))
        mean_long_dev_f.append((pos, np.mean(long_dev_f)))

    for i, (pos, trueFoot) in enumerate(trueFoot_b_arr):
        trueLat, trueLongt = cartesian2latlong(*trueFoot)
        ang_dev_b = []
        lat_dev_b = []
        long_dev_b = []
        (pos, fp_arr) = footpoints_b_arr[i]
        for fp in fp_arr:
            lat, longt = cartesian2latlong(*fp)
            ang_dev_b.append(np.arccos(np.dot(fp, trueFoot)/(norm(fp)*norm(trueFoot))))
            lat_dev_b.append(trueLat - lat)
            long_dev_b.append(trueLongt - longt)
        mean_ang_dev_b.append((pos, np.mean(ang_dev_b)))
        mean_lat_dev_b.append((pos, np.mean(lat_dev_b)))
        mean_long_dev_b.append((pos, np.mean(long_dev_b)))
    
    return mean_ang_dev_f, mean_lat_dev_f, mean_long_dev_f, mean_ang_dev_b, mean_lat_dev_b, mean_long_dev_b
