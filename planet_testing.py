from all_funcs import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib as mat
import numpy as np
from tqdm import tqdm
import os.path
from palettable.wesanderson import Aquatic2_5, Darjeeling2_5


# multiline_3D(10, [0., np.pi/3, 2*np.pi/3], coeffs=uranus)
# plt.show()

def random_footpoints(n, moon, phi):
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
    return footpoints
        


############# TITANIA #############
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
#         ax.plot3D(x, y, z, color=Darjeeling2_5.mpl_colors[3])
#         bar.update()

# print(len(footpoints))
# x, y, z = map(list, zip(*footpoints))
# ax.plot3D(x, y, z, color=Darjeeling2_5.mpl_colors[3])

footpoints = random_footpoints(10, 'Titania', 0)
print(footpoints)
x, y, z  = map(list, zip(*footpoints))
ax.plot3D(x, y, z)

u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:25j]
a = np.cos(u)*np.sin(v)
b = np.sin(u)*np.sin(v)
c = np.cos(v)
# ax.plot_wireframe(a, b, c, color=Aquatic2_5.mpl_colors[0])

plt.show()

