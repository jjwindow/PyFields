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

n_min = 20
n_max = 250
step = 10
# arrs = convergence_test('titania', [17.07, np.pi/2, np.pi], n_min, n_max, step)
_fpaths = ['arange', 'angdevs_f', 'angdevs_b', 'latdev_f', 'latdev_b', 'longdev_f', 'longdev_b']
fpaths = ['Titania/' + path + f'_{n_min}_{n_max}_{step}_phi_pi.npy' for path in _fpaths]
# for arr, path in zip(list(arrs), fpaths):
#     with open(path, 'wb') as file:
#         np.save(file, arr)

with open(fpaths[0], 'rb') as file:
    n_array = np.load(file, allow_pickle=True)
with open(fpaths[1], 'rb') as file:
    angdevs_f = np.load(file, allow_pickle=True)
with open(fpaths[2], 'rb') as file:
    angdevs_b = np.load(file, allow_pickle=True)
with open(fpaths[3], 'rb') as file:
    latdev_f = np.load(file, allow_pickle=True)
with open(fpaths[4], 'rb') as file:
    latdev_b = np.load(file, allow_pickle=True)
with open(fpaths[5], 'rb') as file:
    longdev_f = np.load(file, allow_pickle=True)
with open(fpaths[6], 'rb') as file:
    longdev_b = np.load(file, allow_pickle=True)

# n_array = arrs[0]
# for i in range(1, len(arrs)):
#     plt.plot(n_array, arrs[i], label = _fpaths[i])  

def forward_backward_plots():
    # Make plot of forward and backward footpoint deviations
    # for ang. dev, lat, long. Comment/uncomment as necessary.
    fig, axs = plt.subplots(2, 1, sharex=True)
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xlabel("No. random fieldlines")
    plt.ylabel(r"Ang. deviation squared ($^{\circ}$)")
    # plt.ylabel("Angular deviation (rad)")
    # axs[0].plot(n_array, angdevs_f, label = 'Forward')
    # plt.plot(n_array, angdevs_b, label = 'Backward')

    # GROUP BY LAT-LONG
    # axs[0].set_title("Latitude")
    # axs[0].plot(n_array, [lat**2 for lat in latdev_f], label = "Forward")
    # axs[0].plot(n_array, [lat**2 for lat in latdev_b], label = "Backward")
    # axs[0].legend()
    # axs[1].set_title("Longitude")
    # axs[1].plot(n_array, [long**2 for long in longdev_f], label = "Forward")
    # axs[1].plot(n_array, [long**2 for long in longdev_b], label = "Backward")
    # axs[1].legend()

    # GROUP BY FORWARD-BACKWARD
    axs[0].set_title("Field Into Planet")
    axs[0].plot(n_array, [lat**2 for lat in latdev_f], label = "Latitude")
    axs[0].plot(n_array, [long**2 for long in longdev_f], label = "Longitude")
    axs[0].legend()
    axs[1].set_title("Field Out Of Planet")
    axs[1].plot(n_array, [lat**2 for lat in latdev_b], label = "Latitude")
    axs[1].plot(n_array, [long**2 for long in longdev_b], label = "Longitude")
    axs[1].legend()
    # plt.ylabel("Longitudinal deviation (rad)")
    # plt.plot(n_array, longdev_f, label = "Forward")
    # plt.plot(n_array, longdev_b, label = "Backward")
    plt.show()

def triple_angle_plots():
    # Make plot of ang, lat, long deviations for both forwards
    # and backwards separately.
    fig, axs = plt.subplots(2,1, sharex=True)
    fig.add_subplot(111, frameon=False)
    # hide tick and tick label of the big axis
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.xlabel("Num. Random Fieldlines")
    plt.ylabel("Mean deviation from non-random footpoint (rad)")
    axs[0].plot(n_array, angdevs_f, label="Angle")
    axs[0].plot(n_array, latdev_f, label = "Latitude")
    axs[0].plot(n_array, longdev_f, label = "Longitude")
    axs[0].set_title("Fieldlines into planet")
    axs[0].legend()
    axs[1].plot(n_array, angdevs_b, label="Angle")
    axs[1].plot(n_array, latdev_b, label = "Latitude")
    axs[1].plot(n_array, longdev_b, label = "Longitude")
    axs[1].set_title("Fieldlines out of planet")
    axs[1].legend()

    plt.show()

# triple_angle_plots()
forward_backward_plots()