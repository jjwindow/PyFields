from all_funcs import *
import matplotlib.animation as animation




def trace_partial_orbit(moon, orbit_points, num_fieldlines):
    
    l = len(orbit_points)
    
    footpoints_f_arr = [0 for i in range(l)]
    footpoints_b_arr = [0 for i in range(l)]
    trueFoot_f_arr = [0 for i in range(l)]
    trueFoot_b_arr = [0 for i in range(l)]

    for i, pos in enumerate(orbit_points):
        footpoints_f, footpoints_b, trueFoot_f, trueFoot_b = random_footpoints(num_fieldlines, moon, pos, trueTrace = True)
        footpoints_f_arr[i] = (pos, footpoints_f)
        footpoints_b_arr[i] = (pos, footpoints_b)
        trueFoot_f_arr[i] = (pos, trueFoot_f)
        trueFoot_b_arr[i] = (pos, trueFoot_b)

    return footpoints_f_arr, footpoints_b_arr, trueFoot_f_arr, trueFoot_b_arr


def save_partial_moon_trace(moon, orbit_points, t_s, t_f, num_fieldlines):

    paths = ['footpoints_f', 'footpoints_b', 'trueFoot_f', 'trueFoot_b', 'time']

    # footpoints_f_arr, footpoints_b_arr, trueFoot_f_arr, trueFoot_b_arr = trace_full_orbit(moon, num_orbit_points, num_orbits, num_fieldlines)
    all_footpoints = trace_partial_orbit(moon, orbit_points, num_fieldlines)
    all_footpoints = list(all_footpoints)

    for path, footpoint in zip(paths, all_footpoints):
        fpath = f'{moon}/' + path + f'_partial_{t_s}_{t_f}_{len(orbit_points)}_{num_fieldlines}.npy'
        with open(fpath, 'wb') as file:
            np.save(file, footpoint)


# moon = 'Titania'
# moon = 'Miranda'
# moon = 'Ariel'
def all_moon_footpoints():
    all_moons = ['Miranda', 'Ariel', 'Umbriel', 'Titania', 'Oberon'] 
    ax = plt.subplot(1,1,1)
    for i, moon in enumerate(all_moons):
        if moon == 'Miranda':
            orb = 1.05
        else:
            orb = 1

        with open(f'{moon}/trueFoot_f_40_{orb}_100.npy', 'rb') as file:
            trueFoot_f_arr = np.load(file, allow_pickle=True)
        with open(f'{moon}/trueFoot_b_40_{orb}_100.npy', 'rb') as file:
            trueFoot_b_arr = np.load(file, allow_pickle=True)
        # with open(f'{moon}/footpoints_f_40_{orb}_100.npy', 'rb') as file:
        #     footpoints_f_arr = np.load(file, allow_pickle=True)
        # with open(f'{moon}/footpoints_b_40_{orb}_100.npy', 'rb') as file:
        #     footpoints_b_arr = np.load(file, allow_pickle=True)
        # with open(f'{moon}/time_40_{orb}_100.npy', 'rb') as file:
        #     T_arr = np.load(file, allow_pickle=True)
        for (pos, fp) in trueFoot_f_arr:
            lat, longt = cartesian2latlong(*fp)
            ax.plot(longt + 180, lat, 'x', color=Cavalcanti_5.mpl_colors[i])
        for (pos, fp) in trueFoot_b_arr:
            lat, longt = cartesian2latlong(*fp)
            ax.plot(longt + 180, lat, 'x', color=Cavalcanti_5.mpl_colors[i])
    ax.set_xlabel(r"Longitude ($^{\circ}$)")
    ax.set_ylabel(r"Latitude ($^{\circ}$)")
    ax.set_xlim(360, 0)
    ax.set_ylim(-90, 90)
    ax.legend((mat.lines.Line2D([0,0], [1,1], color = Cavalcanti_5.mpl_colors[0]), mat.lines.Line2D([0,0], [1,1], color = Cavalcanti_5.mpl_colors[1]), mat.lines.Line2D([0,0], [1,1], color = Cavalcanti_5.mpl_colors[2]), mat.lines.Line2D([0,0], [1,1], color = Cavalcanti_5.mpl_colors[3]), mat.lines.Line2D([0,0], [1,1], color = Cavalcanti_5.mpl_colors[4])), ('Miranda', 'Ariel', 'Umbriel', 'Titania', 'Oberon'))
    plt.show()

############ ANIMATE MIRANDA FOOTPOINTS ###########################

fig, ax = plt.subplots()
ax.set_xlabel(r"Longitude ($^{\circ}$)")
ax.set_ylabel(r"Latitude ($^{\circ}$)")
ax.set_xlim(360, 0)
ax.set_ylim(-90, 90)

moon = 'Triton'
ax.set_title(moon)
orb = 2
with open(f'{moon}/trueFoot_f_40_{orb}_100.npy', 'rb') as file:
    trueFoot_f_arr = np.load(file, allow_pickle=True)
with open(f'{moon}/trueFoot_b_40_{orb}_100.npy', 'rb') as file:
    trueFoot_b_arr = np.load(file, allow_pickle=True)
with open(f'{moon}/time_40_{orb}_100.npy', 'rb') as file:
    time = np.load(file, allow_pickle=True)

latf_arr, longf_arr = [], []
latb_arr, longb_arr = [], []
(T,) = moon_selector(moon, 'T')

for (pos, fp_f), (pos, fp_b) in zip(trueFoot_f_arr, trueFoot_b_arr):
    latf, longf = cartesian2latlong(*fp_f)
    latf_arr.append(latf)
    longf_arr.append(longf + 180)
    latb, longb = cartesian2latlong(*fp_b)
    latb_arr.append(latb)
    longb_arr.append(longb + 180)

# time_arr = [(i*1.05)/(len(latf_arr) -1) for i in range(len(latf_arr))]

linef, = ax.plot(longf_arr, latf_arr, '-', label = 'Forward')
lineb, = ax.plot(longb_arr, latb_arr, '-', label = 'Backward')
legend = ax.legend()

def init():
    linef.set_ydata([np.nan] * len(longf_arr))
    lineb.set_ydata([np.nan] * len(longb_arr))
    return linef, lineb,

def animate(i):
    linef.set_xdata(longf_arr[:i])
    linef.set_ydata(latf_arr[:i])
    lineb.set_xdata(longb_arr[:i])
    lineb.set_ydata(latb_arr[:i])
    ax.legend()
    ax.set_title(f"{moon}, t = {round(time[i]/T, 2)}")
    return linef, lineb,

ani = animation.FuncAnimation(fig, animate, frames = len(time), blit = True, init_func = init, interval=400)
ani.save(filename = f'{moon}/{moon}_orbs-{orb}_animation_1.mp4', writer = 'ffmpeg')

####################### WHY THIS THIS ALMOST IDENTICAL TO ABOVE? #####################

# (T,) = moon_selector('Miranda', 'T')
# orbit_points_array, full_time_array = orbit('Miranda', 400, 1, relative=True)
# full_time_array = np.array(full_time_array)/T

def time_filter(t_s, t_f, time_array, orbit_points_array):
    _orb_dict = {t : pos for t, pos in zip(time_array, orbit_points_array)}
    orb_dict = dict(filter(lambda elem: (elem[0] > t_s) and (elem[0] < t_f), _orb_dict.items()))
    # print(list(orb_dict.items())[0], list(orb_dict.items())[-1])
    new_orbit_points_array = list(orb_dict.values())
    new_time_array = list(orb_dict.keys())
    return new_time_array, new_orbit_points_array

# t_s_b = 0.15
# t_f_b = 0.35
# t_s_f = 0.65
# t_f_f = 0.75
# orb = 1.05
# moon = 'Miranda'

# time_f, orbit_points_f = time_filter(t_s_f, t_f_f, full_time_array, orbit_points_array)
# time_b, orbit_points_b = time_filter(t_s_b, t_f_b, full_time_array, orbit_points_array)
# # save_partial_moon_trace('Miranda', orbit_points_f, t_s_f, t_f_f, 1)
# # save_partial_moon_trace('Miranda', orbit_points_b, t_s_b, t_f_b, 1)
# n_points = len(time_b)


# mean_ang_dev_f, mean_lat_dev_f, mean_long_dev_f, mean_ang_dev_b, mean_lat_dev_b, mean_long_dev_b = angular_deviation(footpoints_f_arr, trueFoot_f_arr, footpoints_b_arr, trueFoot_b_arr)
# pos_arr, ang_f = map(list, zip(*mean_ang_dev_f))

# plt.plot(T_arr, ang_f)
# plt.show()

################## PLOT TITANIA FOOTPOINTS FROM SAVED DATA #################

# plt.xlabel(r"Longitude ($^{\circ}$)")
# plt.ylabel(r"Latitude ($^{\circ}$)")
# for (pos, fp) in trueFoot_f_arr:
#     lat, longt = cartesian2latlong(*fp)
#     plt.plot(longt + 180, lat, 'rx')
# for (pos, fp) in trueFoot_b_arr:
#     # print("BACK")
#     lat, longt = cartesian2latlong(*fp)
#     # print((longt, lat))
#     plt.plot(longt + 180, lat, 'bx')
# plt.xlim(360, 0)
# plt.ylim(-90, 90)
# plt.title("Miranda")
# plt.legend((mat.lines.Line2D([0,0], [1,1], color = 'r'), mat.lines.Line2D([0,0], [1,1], color = 'b')), ('Forward', 'Backward'))
# plt.show()

# ####################### PLOT ONE, THEN TWO TRITON ORBITS ############################
# moon = 'Triton'
# orbs = 2

# with open(f'{moon}/trueFoot_f_40_{orbs}_100.npy', 'rb') as file:
#     trueFoot_f_arr = np.load(file, allow_pickle=True)
# with open(f'{moon}/trueFoot_b_40_{orbs}_100.npy', 'rb') as file:
#     trueFoot_b_arr = np.load(file, allow_pickle=True)

# for orb in range(1, orbs+1):
#     fig, ax = plt.subplots()
#     ax.clear()
#     ax.set_xlabel(r"Longitude ($^{\circ}$)")
#     ax.set_ylabel(r"Latitude ($^{\circ}$)")
#     ax.set_xlim(360, 0)
#     ax.set_ylim(-90, 90)
#     ax.set_title(f'{moon}: Orbit {orb}')

#     n = len(trueFoot_f_arr)
#     s = int(np.floor(n*(orb/orbs)))
#     for (pos, fp) in trueFoot_f_arr[:s]:
#         lat, longt = cartesian2latlong(*fp)
#         ax.plot(longt + 180, lat, 'rx')
#     for (pos, fp) in trueFoot_b_arr[:s]:
#         lat, longt = cartesian2latlong(*fp)
#         ax.plot(longt + 180, lat, 'bx')
#     ax.legend((mat.lines.Line2D([0,0], [1,1], color = 'r'), mat.lines.Line2D([0,0], [1,1], color = 'b')), ('Forward', 'Backward'))
#     plt.savefig(f"{moon}_footpoints_40_{orb}_100.png")
#     plt.show()