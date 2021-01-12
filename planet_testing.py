from all_funcs import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mat
import numpy as np
from numpy.linalg import norm
from tqdm import tqdm
import os.path
from palettable.wesanderson import Aquatic2_5, Darjeeling2_5


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

    RETURNS
    ---------------------------------------------------------------------------------
    footpoints - list; list of tuples, where each tuple is (x, y, z) position of a 
                 footpoint of a fieldline calculated from the random coefficients.
    """
    uranian_moons = {'ariel' : [7.469, np.pi/2, phi], 'umbriel' : [10.41, np.pi/2, phi], 
                        'titania' : [17.07, np.pi/2, phi], 'oberon' : [22.83, np.pi/2, phi]}
    # Select coefficients to use 
    moon = moon.lower()
    if moon == 'triton':
        a, g, h = neptune
        a, g_err, h_err = neptune_uncert
        ### START POS NEEDED - MUST BE A POINT ON THE ORBIT. NEED TO FIND GENERIC POINT,
        ### E.G - THE THETA VAL WHERE PHI=0 OR VICE-VERSA
    elif moon in uranian_moons.keys():
        a, g, h = uranus
        a, g_err, h_err = uranus_uncert
        ### SAME CONDITION AS FOR TRITON NEEDED FOR MIRANDA.
        start_pos = uranian_moons[moon]
    else:
        raise ValueError("`moon' arg must be one of the 5 major Uranian moons or 'triton'.")
    
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
        


############# TITANIA #############
# ax = plt.axes(projection = '3d')
# ax.set_xlabel('x')
# ax.set_ylabel('y')
# ax.set_zlabel('z')

# with tqdm(total=50, desc="FOOTPOINTS") as bar:
#     footpoints = []
#     for phi in np.linspace(0, 2*np.pi, 50):
#         start_pos = [17.188, (np.pi/2 - 0.00593), phi]
#         x, y, z = field_trace(start_pos, uranus, 0.005, 200000)
#         # point = (x[-1], y[-1], z[-1])
#         # footpoints.append(point)
#         ax.plot3D(x, y, z, color=Darjeeling2_5.mpl_colors[3])
#         bar.update()

# print(len(footpoints))
# x, y, z = map(list, zip(*footpoints))
# ax.plot3D(x, y, z, color=Darjeeling2_5.mpl_colors[3])

# u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:25j]
# a = np.cos(u)*np.sin(v)
# b = np.sin(u)*np.sin(v)
# c = np.cos(v)
# # ax.plot_wireframe(a, b, c, color=Aquatic2_5.mpl_colors[0])

# plt.show()

###### Plotting range of footpoints for a single position on lat-long plot ######
phi = 0
moon = 'Titania'
# footpoints, trueFoot = random_footpoints(100, moon, phi, True)
fpath = 'Titania_phi-0_n-100.npy'
fpathTrue = 'Titania_phi-0_true.npy'
# with open(fpath, 'wb') as file:
#     np.save(file, footpoints)
# with open(fpathTrue, 'wb') as file:
#     np.save(file, np.asarray(trueFoot))

with open(fpath, 'rb') as file:
    footpoints = np.load(file, allow_pickle=True)
with open(fpathTrue, 'rb') as file:
    trueFoot = np.load(file, allow_pickle=True)
    
x, y, z  = map(list, zip(*footpoints))
lat, longt = cartesian2latlong(x, y, z)
trueLat, trueLongt = cartesian2latlong(*trueFoot)


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
ang_dev = []
lat_dev = []
long_dev = []
for p in footpoints:
    p = np.asarray(p)
    lat, longt = cartesian2latlong(*p)
    ang_dev.append(np.arccos(np.dot(p, trueFoot)/(norm(p)*norm(trueFoot))))
    lat_dev.append(trueLat - lat)
    long_dev.append(trueLongt - longt)

mean_ang_dev = np.mean(ang_dev)
mean_lat_dev = np.mean(lat_dev)
mean_long_dev = np.mean(long_dev)
print("Mean Ang: ", mean_ang_dev*180/np.pi)
print("Mean lat: ", mean_lat_dev)
print("Mean long: ", mean_long_dev)

# plt.hist(ang_dev, bins='auto', edgecolor='k',)
# plt.xlabel(r"Angular Deviation of Footpoints ($^\circ$)")
# plt.ylabel("Frequency Density")
# plt.axvline(mean_ang_dev, color = 'k', linestyle='dashed', label="Arithmetic Mean")
# plt.legend()
# plt.show()


###### Histograms ######
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


###### Plotting an entire orbit on a lat-long plot ######
# with tqdm(total=50, desc="FOOTPOINTS") as bar:
#     latitudes = []
#     longitudes = []
#     for phi in np.linspace(0, 2*np.pi, 100):
#         start_pos = [17.07, np.pi/2, phi]
#         x, y, z = field_trace(start_pos, uranus, 0.005, 200000)
#         latitude, longitude = cartesian2latlong(x[-1], y[-1], z[-1])
#         latitudes.append(latitude)
#         longitudes.append(longitudes)
#         plt.plot(longitudes, latitudes, 'x')
#         bar.update()

# plt.xlabel(r'Longitude ($^\circ$)')
# plt.ylabel(r'Latidude ($^\circ$)')
# plt.show()

def orbit(planet, radius, period_moon, period_plan, incl, num, num_orbits):      #num_orbits is how many sidereal orbits #num gives num of points in one sidereal orbit
    omega_moon = (2*np.pi)/period_moon
    omega_plan = (2*np.pi)/period_plan
    t_step = period_moon/num
    n = num*num_orbits

    footpoints = []

    for n in np.linspace(0, n, n+1):
        pos = [radius, np.pi/2 - incl*np.sin(omega_moon*n*t_step), (omega_moon*n*t_step)-(omega_plan*n*t_step)]
        x, y, z = field_trace(pos, planet, 0.005, 200000)
        point = (x[-1], y[-1], z[-1])
        footpoints.append(points)
    
    return footpoints

