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
    if moon == 'Triton':
        a, g, h = neptune
        a, g_err, h_err = neptune_uncert
    else:
        a, g, h = uranus
        a, g_err, h_err = uranus_uncert
        if moon == 'Titania':
            start_pos = [17.188, np.pi/2, phi]
        else: 
            pass

    footpoints = [0. for _ in range(n)]
    with tqdm(total=n, desc=":)") as bar:
        for k in range(n):
            g_new = np.zeros((3,3))
            h_new = np.zeros((3,3))

            for i in range(3):
                for j in range(3):
                    if g[i][j] == 0.:
                        pass
                    else:
                        r = (np.random.random()-0.5)*2
                        g_new[i][j] = g[i][j] + g_err[i][j]*r
                        h_new[i][j] = h[i][j] + h_err[i][j]*r
            
            coeffs = (a, g_new, h_new)
            x, y, z = field_trace(start_pos, coeffs, 0.005, 200000)
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

footpoints = random_footpoints(100, 'Titania', 0)
x, y, z  = map(list, zip(*footpoints))
ax.plot3D(x, y, z)

u, v = np.mgrid[0:2*np.pi:50j, 0:np.pi:25j]
a = np.cos(u)*np.sin(v)
b = np.sin(u)*np.sin(v)
c = np.cos(v)
# ax.plot_wireframe(a, b, c, color=Aquatic2_5.mpl_colors[0])

plt.show()

