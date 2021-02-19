from all_funcs import *
from matplotlib import animation

all_moons = ['Miranda', 'Ariel', 'Umbriel', 'Titania', 'Oberon', 'Triton'] 

def plot_ang_dev_timeseries(moon, trueFoot_f, trueFoot_b, footpoints_f, footpoints_b, time):
    angs_array = angular_deviation(footpoints_f, trueFoot_f, footpoints_b, trueFoot_b)
    fignames = [f"{moon} Forwards Angular Deviation", f"{moon} Forwards Latitude Deviation", f"{moon} Forwards Longitude Deviation", f"{moon} Backwards Angular Deviation", "Backwards Latitude Deviation", "Backwards Longitude Deviation"]
    T_m, T_p = moon_selector(moon, 'T', 'parent_day')
    T_rel = abs(T_m*T_p/(T_m-T_p))
    time = time/T_rel
    ax = plt.subplot(1,1,1)
    n = int(len(fignames)/2)
    for i, (f_array, name) in enumerate(zip(angs_array, fignames)):
        ax.clear()
        if i >= n:
            break
        f_array = [ang for (pos, ang) in f_array]
        ax.plot(time, f_array, label = name)
        b_array = angs_array[i+3]
        b_array = [ang for (pos, ang) in b_array]
        ax.plot(time, b_array, label = fignames[i+3])

        ax.legend()
        plt.show()

def ang_devs_timeseries(moon_arr):

    for i, moon in enumerate(moon_arr):

        if moon == 'Miranda':
            orb = 1.05
        elif moon == 'Triton':
            orb = 2
        else:
            orb = 1

        with open(f'{moon}/trueFoot_f_40_{orb}_100.npy', 'rb') as file:
            trueFoot_f_arr = np.load(file, allow_pickle=True)
        with open(f'{moon}/trueFoot_b_40_{orb}_100.npy', 'rb') as file:
            trueFoot_b_arr = np.load(file, allow_pickle=True)
        with open(f'{moon}/footpoints_f_40_{orb}_100.npy', 'rb') as file:
            footpoints_f_arr = np.load(file, allow_pickle=True)
        with open(f'{moon}/footpoints_b_40_{orb}_100.npy', 'rb') as file:
            footpoints_b_arr = np.load(file, allow_pickle=True)
        with open(f'{moon}/time_40_{orb}_100.npy', 'rb') as file:
            time = np.load(file, allow_pickle=True)

        plot_ang_dev_timeseries(moon, trueFoot_f_arr, trueFoot_b_arr, footpoints_f_arr, footpoints_b_arr, time)

# ang_devs_timeseries(all_moons)

##################### DATA RETRIEVAL ##########################

def get_coarse_moon_footpoints(moon, n_o_p, n_o):
    """
    Get specific data for moon with n_o_p number points per orbit and n_o 
    number of orbits.
    """
    with open(f'{moon}/trueFoot_f_{n_o_p}_{n_o}_100.npy', 'rb') as file:
        trueFoot_f_arr = np.load(file, allow_pickle=True)
    with open(f'{moon}/trueFoot_b_{n_o_p}_{n_o}_100.npy', 'rb') as file:
        trueFoot_b_arr = np.load(file, allow_pickle=True)
    with open(f'{moon}/time_{n_o_p}_{n_o}_100.npy', 'rb') as file:
        time = np.load(file, allow_pickle=True)
    return trueFoot_f_arr, trueFoot_b_arr, time

def get_final_moon_footpoints(moon):
    """
    Gets forward, backward and time arrays from final (n=500) runs.
    """
    if moon == 'Triton':
        num_orbits = 2
    else:
        num_orbits = 1
    with open(f"Finals/{moon}/trueFoot_f_500_{num_orbits}.npy", 'rb') as file:
        trueFoot_f_arr = np.load(file, allow_pickle=True)
    with open(f"Finals/{moon}/trueFoot_b_500_{num_orbits}.npy", 'rb') as file:
        trueFoot_b_arr = np.load(file, allow_pickle=True)
    with open(f'Finals/{moon}/time_500_{num_orbits}.npy', 'rb') as file:
        time = np.load(file, allow_pickle=True)
    return trueFoot_f_arr, trueFoot_b_arr, time

################# ANIMATIONS #######################

def animate_footpoints(moon, footpoints, delay):
    """
    footpoints = (trueFoot_f_arr, trueFoot_b_arr, time)
    delay (ms) - interval passed to animate
    """
    (trueFoot_f_arr, trueFoot_b_arr, time) = footpoints
    P = max(time)
    # Setup figure axes
    fig, ax = plt.subplots()
    ax.set_xlabel(r"Longitude ($^{\circ}$)")
    ax.set_ylabel(r"Latitude ($^{\circ}$)")
    ax.set_xlim(360, 0)
    ax.set_ylim(-90, 90)

    latf_arr, longf_arr = [], []
    latb_arr, longb_arr = [], []

    # Calc latitudes + longitudes
    for (pos, fp_f), (pos, fp_b) in zip(trueFoot_f_arr, trueFoot_b_arr):
        latf, longf = cartesian2latlong(*fp_f)
        latf_arr.append(latf)
        longf_arr.append(longf + 180)
        latb, longb = cartesian2latlong(*fp_b)
        latb_arr.append(latb)
        longb_arr.append(longb + 180)

    linef, = ax.plot(longf_arr, latf_arr, '-', label = 'Forward')
    lineb, = ax.plot(longb_arr, latb_arr, '-', label = 'Backward')
    legend = ax.legend()

    def init():
        # Flush figures to begin
        linef.set_ydata([np.nan] * len(longf_arr))
        lineb.set_ydata([np.nan] * len(longb_arr))
        return linef, lineb,

    def animate(i):
        linef.set_xdata(longf_arr[:i])
        linef.set_ydata(latf_arr[:i])
        lineb.set_xdata(longb_arr[:i])
        lineb.set_ydata(latb_arr[:i])
        ax.legend()
        ax.set_title(f"{moon}, t = {round(time[i]/(P), 2)}" + r" x $T_{rel}$")
        return linef, lineb,

    ani = animation.FuncAnimation(fig, animate, frames = len(time), blit = True, init_func = init, interval=delay)
    return ani

def animate_all_moons():
    print("Animating Moons...")
    for moon in all_moons:
        print(f"...{moon}")
        n_o = (moon == 'Triton')*2 + (moon != 'Triton')
        fpath = f'Finals/{moon}/Animation_{moon}_500_{n_o}.mp4'
        footpoints = get_final_moon_footpoints(moon)
        ani = animate_footpoints(moon, footpoints, 50)
        ani.save(filename = fpath, writer = 'ffmpeg')
    print("Done")

# animate_all_moons()

######## COORDINATE TIMESERIES ########

def lat_long_timeseries(moon):
    
    trueFoot_f_arr, trueFoot_b_arr, time = get_final_moon_footpoints(moon)
    T_rel = max(abs(time))
    time /= T_rel

    latlongs_f = [cartesian2latlong(*foot) for (pos, foot) in trueFoot_f_arr]
    lat_f, long_f = map(np.asarray, zip(*latlongs_f))
    long_f += 180

    latlongs_b = [cartesian2latlong(*foot) for (pos, foot) in trueFoot_b_arr]
    lat_b, long_b = map(np.asarray, zip(*latlongs_b))
    long_b += 180

    fig, axs = plt.subplots(2, 1, sharex=True)
    axs[0].clear()
    axs[1].clear()
    axs[0].set_ylim(0, 360)
    axs[0].plot(time, long_f, label='Forwards')
    axs[0].plot(time, long_b, label = 'Backwards')
    axs[0].set_ylabel(r"Longitude ($^{\circ}$)")
    
    axs[1].set_ylim(-90, 90)
    axs[1].plot(time, lat_f, label='Forwards')
    axs[1].plot(time, lat_b, label = 'Backwards')
    axs[1].set_ylabel(r"Latitude ($^{\circ}$)")
    axs[1].set_xlabel(r"Time ($T_{rel}$)")
    axs[0].set_title(f"{moon}")

    axs[0].legend()
    axs[1].legend()
    
    plt.savefig(f"Finals/{moon}/timeseries_{moon}_lat_long.png")
    print(f"{moon} Figure saved.")

[lat_long_timeseries(moon) for moon in all_moons]

def footpoint_velocities(moon):

    trueFoot_f_arr, trueFoot_b_arr, time = get_final_moon_footpoints(moon)
    T_rel = max(abs(time))
    time /= T_rel



# fig, axs = plt.subplots(2, 1, sharex=True)
# n = int(1.5*len(time)/5)
# axs[0].clear()
# axs[1].clear()
# axs[0].plot(time[:n]/(0.2*T_m*n_o), latf_arr[:n], 'b-', label = 'Forwards')
# axs[0].plot(time[:n]/(0.2*T_m*n_o), latb_arr[:n], 'r-', label = 'Backwards')
# axs[0].set_ylabel(r"Latitude ($^{\circ}$)")
# axs[1].plot(time[:n]/(0.2*T_m*n_o), longf_arr[:n], 'b-', label = 'Forwards')
# axs[1].plot(time[:n]/(0.2*T_m*n_o), longb_arr[:n], 'r-', label = 'Backwards')
# axs[1].set_ylabel(r"Longitude ($^{\circ}$)")
# axs[1].set_xlabel(r"Time /$T_{rel}$")
# axs[0].legend()
# axs[1].legend()

# # plt.show()

######## LOWES SPECTRA & FIELD RATIO DECAY ###########################

# r = np.linspace(1, 25, 1000)
# ratio = r**-1
# quad = r**-4
# dip = r**-3

# uranus_Lowes = np.array([(i+1)*sum(_g**2 for _g in g) for i, g in enumerate(g_U)]) + np.array([(i+1)*sum(_h**2 for _h in h) for i, h in enumerate(h_U)])
# uranus_Lowes = uranus_Lowes[1:]
# uranus_Lowes /= uranus_Lowes[0]

# neptune_Lowes = np.array([(i+1)*sum(_g**2 for _g in g) for i, g in enumerate(g_N)]) + np.array([(i+1)*sum(_h**2 for _h in h) for i, h in enumerate(h_N)])
# neptune_Lowes = neptune_Lowes[1:]
# neptune_Lowes /= neptune_Lowes[0]

# ratio_U = uranus_Lowes[1]/uranus_Lowes[0]
# ratio_N = neptune_Lowes[1]/neptune_Lowes[0]




# # plt.plot(range(len(uranus_Lowes)), uranus_Lowes)
# # plt.plot(range(len(neptune_Lowes)), neptune_Lowes)

# plt.plot(r, ratio_U*ratio, label = 'Uranus')
# plt.plot(r, ratio_N*ratio, label = 'Neptune')
# for moon in all_moons:
#     a, = moon_selector(moon, 'a')
#     plt.axvline(a, linestyle = '-.', color = 'k')
#     plt.text(a+0.5, 1, f'{moon}', rotation=90)
# plt.legend()
# plt.show()

