from all_funcs import *
from all_funcs import _B_rad_


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

def time_filter(t_s, t_f, time_array, orbit_points_array):
    _orb_dict = {t : pos for t, pos in zip(time_array, orbit_points_array)}
    orb_dict = dict(filter(lambda elem: (elem[0] > t_s) and (elem[0] < t_f), _orb_dict.items()))
    # print(list(orb_dict.items())[0], list(orb_dict.items())[-1])
    new_orbit_points_array = list(orb_dict.values())
    new_time_array = list(orb_dict.keys())
    return new_time_array, new_orbit_points_array

def surface_radial_field(planet):

    coeffs = (planet == 'Uranus')*uranus + (planet == 'Neptune')*neptune
    # rv, thv, phv = np.meshgrid(np.ones(500), np.linspace(0, np.pi, 500), np.linspace(0, 2*np.pi, 1000))
    # points = [(r, th, ph) for r, (th, ph) in zip(np.ones(10), zip(np.linspace(0, np.pi, 10), np.linspace(0, 2*np.pi, 10)))]
    points = [[(th, ph) for ph in np.linspace(0, 2*np.pi, 500)] for th in np.linspace(0, np.pi, 500)]
    # B_r_mat = _B_rad_(rv, thv, phv, *coeffs)
    B_r_mat = [[_B_rad_(1, *p, *coeffs) for p in row] for row in points]
    plt.imshow(B_r_mat, cmap='magma', extent = [360, 0, -90, 90])
    plt.colorbar()
    plt.show()

# surface_radial_field('Uranus')

def surface_radial_field_contour(planet):

    coeffs = (planet == 'Uranus')*uranus + (planet == 'Neptune')*neptune
    # rv, thv, phv = np.meshgrid(np.ones(500), np.linspace(0, np.pi, 500), np.linspace(0, 2*np.pi, 1000))
    # points = [(r, th, ph) for r, (th, ph) in zip(np.ones(10), zip(np.linspace(0, np.pi, 10), np.linspace(0, 2*np.pi, 10)))]
    points = [[(th, ph) for ph in np.linspace(0, 2*np.pi, 500)] for th in np.linspace(0, np.pi, 500)]
    # B_r_mat = _B_rad_(rv, thv, phv, *coeffs)
    B_r_mat = np.asarray([np.asarray([_B_rad_(1, *p, *coeffs) for p in row]) for row in points])
    print(np.amin(B_r_mat))
    B_r_mat += abs(np.amin(B_r_mat))
    print(np.amin(B_r_mat))
    fig, ax = plt.subplots()
    cs = ax.contourf(np.linspace(360, 0, 500), np.linspace(-90, 90, 500), B_r_mat)
    ax.invert_xaxis()
    # plt.imshow(B_r_mat, cmap='magma', extent = [360, 0, -90, 90])
    plt.colorbar(cs)
    plt.show()

# surface_radial_field('Neptune')