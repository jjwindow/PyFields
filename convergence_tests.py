from all_funcs import *


def convergence_test(moon, pos, n_min, n_max, step):
    n_array = np.arange(n_min, n_max+1, step)
    ang_dev_arr_f = []
    ang_dev_arr_b = []
    lat_dev_arr_f = []
    lat_dev_arr_b = []
    long_dev_arr_f = []
    long_dev_arr_b = []

    for n in n_array:
        footpoints_f, footpoints_b, trueFoot_f, trueFoot_b = random_footpoints(n, moon, pos, True)
        mean_ang_dev_f, mean_lat_dev_f, mean_long_dev_f, mean_ang_dev_b, mean_lat_dev_b, mean_long_dev_b = angular_deviation([(pos, footpoints_f)], [(pos, trueFoot_f)], [(pos, footpoints_b)], [(pos, trueFoot_b)])
        (pos, ang_dev_f) = mean_ang_dev_f[0]
        ang_dev_arr_f.append(ang_dev_f)
        (pos, ang_dev_b) = mean_ang_dev_b[0]
        ang_dev_arr_b.append(ang_dev_b)
        (pos, lat_dev_f) = mean_lat_dev_f[0]
        lat_dev_arr_f.append(lat_dev_f)
        (pos, lat_dev_b) = mean_lat_dev_b[0]
        lat_dev_arr_b.append(lat_dev_b)
        (pos, long_dev_f) = mean_long_dev_f[0]
        long_dev_arr_f.append(long_dev_f)
        (pos, long_dev_b) = mean_long_dev_b[0]
        long_dev_arr_b.append(long_dev_b)
    
    return n_array, ang_dev_arr_f, ang_dev_arr_b, lat_dev_arr_f, lat_dev_arr_b, long_dev_arr_f, long_dev_arr_b

n_min = 5
n_max = 500
step = 10
arrs = convergence_test('titania', [17.07, np.pi/2, np.pi], 5, 500, 10)
fpaths = ['arange', 'angdevs_f', 'angdevs_b', 'latdev_f', 'latdev_b', 'longdev_f', 'longdev_b']
fpaths = ['Titania/' + path + f'_{n_min}_{n_max}_{step}_phi_pi.npy' for path in fpaths]
for arr, path in zip(list(arrs), fpaths):
    with open(path, 'wb') as file:
        np.save(file, n_array)
        
plt.plot(n_array, ang_dev_arr_f, 'bx')
plt.show()
