from all_funcs import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mat
import numpy as np
from numpy.linalg import norm
from tqdm import tqdm
import os.path
from palettable.wesanderson import Aquatic2_5, Cavalcanti_5


# multiline_3D(10, [0., np.pi/3, 2*np.pi/3], coeffs=uranus)
# plt.show()

def random_footpoints(n, moon, phi, trueTrace = False):
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

    (R, coeffs, uncert) = moon_selector(moon, 'R', 'coeffs', 'uncert')
    start_pos = [R, np.pi/2, phi]
    (a, g, h) = coeffs
    
    # Trace the accepted fieldline if desired
    if trueTrace:
        x, y, z = field_trace(start_pos, (a, g, h), 0.005, 200000)
        trueFoot = (x[-1], y[-1], z[-1])
    # initialise footpoints array
    footpoints = [0. for _ in range(n)]
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
            x, y, z = field_trace(start_pos, coeffs, 0.005, 200000)
            # Take fieldline footpoint
            footpoints[k] = (x[-1], y[-1], z[-1])
            bar.update()

    if trueTrace:
        return footpoints, trueFoot
    else:
        return footpoints
        

###### Plotting range of footpoints for a single position on lat-long plot ######
# phi = 0
# moon = 'Titania'
# footpoints, trueFoot = random_footpoints(100, moon, phi, True)
# fpath = 'Titania_phi-0_n-100.npy'
# fpathTrue = 'Titania_phi-0_true.npy'
# with open(fpath, 'wb') as file:
#     np.save(file, footpoints)
# with open(fpathTrue, 'wb') as file:
#     np.save(file, np.asarray(trueFoot))

# with open(fpath, 'rb') as file:
#     footpoints = np.load(file, allow_pickle=True)
# with open(fpathTrue, 'rb') as file:
#     trueFoot = np.load(file, allow_pickle=True)
    
# x, y, z  = map(list, zip(*footpoints))
# lat, longt = cartesian2latlong(x, y, z)
# trueLat, trueLongt = cartesian2latlong(*trueFoot)


def makeThisAPlottingFunc():
    """
    Make this a general plotting func later.
    """
    plt.plot(trueLongt, trueLat, 'ro', label = r"Accepted $g_n^m,~h_n^m$")
    plt.plot(longt, lat, 'x', label = r"Random $g_n^m,~h_n^m$")
    plt.annotate(f"{moon}, phi = {phi}", (0.7, 0.05), xycoords = 'axes fraction')
    plt.xlabel(r'Longitude ($^\circ$)')
    plt.ylabel(r'Latidude ($^\circ$)')
    plt.legend()
    plt.show()

###### Calculating mean angular error ######
# ang_dev = []
# lat_dev = []
# long_dev = []
# for p in footpoints:
#     p = np.asarray(p)
#     lat, longt = cartesian2latlong(*p)
#     ang_dev.append(np.arccos(np.dot(p, trueFoot)/(norm(p)*norm(trueFoot))))
#     lat_dev.append(trueLat - lat)
#     long_dev.append(trueLongt - longt)

# mean_ang_dev = np.mean(ang_dev)
# mean_lat_dev = np.mean(lat_dev)
# mean_long_dev = np.mean(long_dev)
# print("Mean Ang: ", mean_ang_dev*180/np.pi)
# print("Mean lat: ", mean_lat_dev)
# print("Mean long: ", mean_long_dev)


###### Histograms ######
def histograms_dep():
    """
    Histograms of angular deviation due to uncertainty of harmonic coefficients.
    """
    lat_devs = []
    longt_devs = []
    latitudes = []
    longitudes = []
    for fp in footpoints:
        x, y, z = fp
        latitude, longitude = cartesian2latlong(x, y, z)
        latitudes.append(latitude)
        longitudes.append(longitude)
        lat_devs.append(trueLat - latitude)
        longt_devs.append(trueLongt - longitude)

    # fig, ax1 = plt.subplots(3, 1, sharex = True)
    fig = plt.figure()
    title_ax = fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    title_ax.set_ylabel("Frequency Density")
    title_ax.set_xlabel(r"Deviation from Accepted Footpoint ($^\circ$)")
    ax1 = fig.add_subplot(3,1,1)
    ax1.hist(lat_devs, bins='auto', color='b', edgecolor='k', label="Latitude")  #latitude deviations histogram
    ax1.axvline(mean_lat_dev, color='k', linestyle='dashed', linewidth=1, label = f"Mean: {round(mean_lat_dev, 3)}")
    ax1.legend()

    ax2 = fig.add_subplot(3,1,2)
    ax2.hist(longt_devs, bins='auto', color='c', edgecolor='k', label = "Longitude")   #longitude deviations histogram
    ax2.axvline(mean_long_dev, color='k', linestyle='dashed', linewidth=1, label = f"Mean: {round(mean_long_dev, 3)}")
    ax2.legend()

    ax3 = fig.add_subplot(313)
    ax3.hist([ang*180/np.pi for ang in ang_dev], bins='auto', edgecolor='k', label = 'Absolute Angle')
    ax3.axvline(mean_ang_dev*180/np.pi, color = 'k', linestyle='dashed', label=f"Mean: {round(mean_ang_dev*180/np.pi, 3)}")
    ax3.legend()
    plt.show()

# fig, ax2 = plt.subplots(2, 1)
# ax2[0].hist(latitudes, bins='auto')   #latitudes histogram
# ax2[0].axvline(trueLat, color='k', linestyle='dashed', linewidth=1)

# ax2[1].hist(longitudes, bins='auto')   #longitudes histogram
# ax2[1].axvline(trueLongt, color='k', linestyle='dashed', linewidth=1)


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

    for i in range(n+1):
        t = i * t_step      # elapsed time
        # angular argument of satellite in the plane of its orbit, more correctly called the 'argument of latitude'.
        phi_moon_orbit = omega_moon * t     
        # from Adam's eqns:
        theta = np.arccos(np.cos(phi_moon_orbit)*np.sin(np.pi-incl))    
        phi_moon_eq = np.arctan2(-1*np.sin(phi_moon_orbit), np.cos(phi_moon_orbit)*np.cos(np.pi - incl))
        # phi_moon_eq is latitude coordinate in equatorial plane.
        if phi_moon_eq < 0:
            # handles negative arctan2 output
            print(phi_moon_eq)
            phi_moon_eq += 2*np.pi
        if relative:
            # changes to planet rest frame
            phi = phi_moon_eq - omega_plan * t
        else:
            phi = phi_moon_eq 
        # append point to list
        pos = [R, theta, phi]
        orbital_points[i] = pos
    return np.array(orbital_points)

############# ORBIT TESTING #############

# set up 3d axes
ax = plt.axes(projection = '3d')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

# with tqdm(total=50, desc="FOOTPOINTS") as bar:
#     footpoints = []
#     for phi in np.linspace(0, 2*np.pi, 50):
#         start_pos = [17.188, (np.pi/2 - 0.00593), phi]
#         x, y, z = field_trace(start_pos, uranus, 0.005, 200000)
#         # point = (x[-1], y[-1], z[-1])
#         # footpoints.append(point)
#         ax.plot3D(x, y, z, color=Cavalcanti_5.mpl_colors[3])
#         bar.update()

# print(len(footpoints))
# x, y, z = map(list, zip(*footpoints))

def plot_orbits(moons_list, relative = False):
    """
    Plots all orbital paths for moons in a 'moons_list', in sidereal or planet
    rest frame, on 3d axes. Also plots planet for scale.
    PARAMS
    -----------------------------------------------------------------------------
    moons_list  -   array or list; contains only elements of type str, which must
                    be one of the Uranian or Neptunian moons.
    relative    -   bool; if false, orbits plotted in sidereal rest frame. Otherwise,
                    plotted in planet rest frame.
    """
    # initialise lists
    x_ptp_arr = []
    y_ptp_arr = []
    z_ptp_arr = []

    # plot each moon in list
    for i, moon in enumerate(moons_list):
        orbital_points = orbit(moon, 50, 1)     # retrieve orbital path
        x, y, z = spherical2cartesian(orbital_points)       # convert to Cartesian
        ax.plot3D(x, y, z, color=Cavalcanti_5.mpl_colors[i], label = moon)
        # save peak-to-peak width of orbital path in each co-ord.
        x_ptp_arr.append(np.ptp(x))
        y_ptp_arr.append(np.ptp(y))
        z_ptp_arr.append(np.ptp(z))

    # plot planet.
    u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:25j]
    a = np.cos(u)*np.sin(v)
    b = np.sin(u)*np.sin(v)
    c = np.cos(v)
    # find maximum bound for each coordinate
    x_len = max(x_ptp_arr)
    y_len = max(y_ptp_arr)
    # maximum z-bound either set by orbit or by planet, must compare both
    z_len = max(max(z_ptp_arr), np.ptp(c))
    # set aspect ratio by largest path in each dimension -> no squished paths or planets
    ax.set_box_aspect((x_len, y_len, z_len))
    ax.plot_wireframe(a, b, c, color=Aquatic2_5.mpl_colors[0])
    plt.legend()

    plt.show()

# Plotting the different planetary systems
uranus_moons = ['Miranda', 'Ariel', 'Umbriel', 'Titania', 'Oberon']
# plot_orbits(uranus_moons)
plot_orbits(['Triton'], True)