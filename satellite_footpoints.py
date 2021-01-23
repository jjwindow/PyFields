from all_funcs import *


def trace_full_orbit(moon, num_orbit_points, num_orbits, num_fieldlines):
    
    orbital_points_arr = orbit(moon, num_orbit_points, num_orbits, relative = True)
    l = len(orbital_points_arr)
    
    footpoints_f_arr = [0 for i in range(l)]
    footpoints_b_arr = [0 for i in range(l)]
    trueFoot_f_arr = [0 for i in range(l)]
    trueFoot_b_arr = [0 for i in range(l)]

    for i, pos in enumerate(orbital_points_arr):
        footpoints_f, footpoints_b, trueFoot_f, trueFoot_b = random_footpoints(num_fieldlines, moon, pos, trueTrace = True)
        footpoints_f_arr[i] = (pos, footpoints_f)
        footpoints_b_arr[i] = (pos, footpoints_b)
        trueFoot_f_arr[i] = (pos, trueFoot_f)
        trueFoot_b_arr[i] = (pos, trueFoot_b)

    return footpoints_f_arr, footpoints_b_arr, trueFoot_f_arr, trueFoot_b_arr

def save_moon_trace(moon, num_orbit_points, num_orbits, num_fieldlines):

    paths = ['footpoints_f', 'footpoints_b', 'trueFoot_f', 'trueFoot_b']

    # footpoints_f_arr, footpoints_b_arr, trueFoot_f_arr, trueFoot_b_arr = trace_full_orbit(moon, num_orbit_points, num_orbits, num_fieldlines)
    all_footpoints = trace_full_orbit(moon, num_orbit_points, num_orbits, num_fieldlines)
    all_footpoints = list(all_footpoints)

    for path, footpoint in zip(paths, all_footpoints):
        fpath = f'{moon}/' + path + f'_{num_orbit_points}_{num_orbits}_{num_fieldlines}.npy'
        with open(fpath, 'wb') as file:
            np.save(file, footpoint)

# save_moon_trace('Titania', 10, 1, 10)

with open('Titania/trueFoot_f_10_1_10.npy', 'rb') as file:
    trueFoot_f_arr = np.load(file, allow_pickle=True)
with open('Titania/trueFoot_b_10_1_10.npy', 'rb') as file:
    trueFoot_b_arr = np.load(file, allow_pickle=True)

all_moons = ['Miranda', 'Ariel', 'Umbriel', 'Titania', 'Oberon', 'Triton']
plt.xlabel(r"Longitude ($^{\circ}$)")
plt.ylabel(r"Latitude ($^{\circ}$)")
for (pos, fp) in trueFoot_f_arr:
    lat, longt = cartesian2latlong(*fp)
    plt.plot(longt, lat, 'rx')
for (pos, fp) in trueFoot_b_arr:
    print("BACK")
    lat, longt = cartesian2latlong(*fp)
    print((longt, lat))
    plt.plot(longt, lat, 'bx')
plt.legend((mat.lines.Line2D([0,0], [1,1], color = 'r'), mat.lines.Line2D([0,0], [1,1], color = 'b')), ('Forward', 'Backward'))
plt.show()
