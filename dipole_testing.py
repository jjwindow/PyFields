"""
PyFields dipole_testing.py

28.10.2020

Testing model using a dipole
"""

from all_funcs import *
import matplotlib.pyplot as plt
import matplotlib as mat
import numpy as np
from tqdm import tqdm
import os.path

########################## DIPOLE TEST ###############################

# multiline_plot(50)
# params = {
#  'axes.labelsize': 14,
#    'font.size': 14,
#    'legend.fontsize': 14,
#    'xtick.labelsize': 12,
#    'ytick.labelsize': 12,
#    'figure.figsize': [8,6]
#    }
# plt.rcParams.update(params)
# plt.legend((mat.lines.Line2D([0,0], [1,1], color = 'r'),), ('Traced Dipole, ds = 0.01',))
# plt.xlabel("Distance in x (Scaled by Planetary Radius)")
# plt.ylabel("Distance in y (Scaled by Planetary Radius)")
# plt.show()


########################## ERROR CALCULATION ##########################

def dipole_error(num, th_min, th_max, ds, max_iter):
    th_values = np.linspace(th_min, th_max, num, endpoint=False)
    th_returns = []
    for th in th_values:
        if (th < 0.0):
            th_return = -np.pi + abs(th)
            th_returns.append(th_return)
        else:
            th_return = np.pi - th
            th_returns.append(th_return)
    
    th_returns = np.array(th_returns)

    deltas = []
    lengths = []
    th_finals = []
    with tqdm(total=len(th_values), desc="Tracing Fields...") as bar:
        for i, th in enumerate(th_values):
            start_pos = [1., th, 0.]
            field = field_trace(start_pos, dipole, ds, max_iter, axes=None)
            if field is not None:
                (p_arr, B_arr) = field
                th_final = p_arr[-1][1]
                th_finals.append(th_final)
                lengths.append(len(p_arr))
                deltas.append(abs(th_final-th_returns[i]))
            else:
                th_values[i] = np.nan
                th_finals.append(np.nan)
                lengths.append(np.nan)
                deltas.append(np.nan)
            bar.update()
       
    deltas = np.array(deltas)
    lengths = np.array(lengths)

    return th_values, deltas, lengths

def multi_step_size(num, th_min, th_max, stepsizes):
    for ds in stepsizes:
        field_lines = multiline_plot(num, th_min, th_max, ds=ds, maxits= int(ds*1e7), plot=False)
        th_values, deltas, lengths = dipole_error(num, th_min, th_max, ds, int(ds*1e7))
        fpath_field = f'/Testing/Dipole/Fieldlines/Dipole_fieldlines_ds_{ds}.npy'
        fpath_errors = f'/Testing/Dipole/Errors/Dipole_errors_ds_{ds}.npy'
        with open(fpath_field, 'wb') as file:
            np.save(file, field_lines)
        with open(fpath_errors, 'wb') as file:
            np.save(file, [th_values, deltas, lengths])

stepsizes=[0.1, 0.01]
multi_step_size(50, -np.pi/2, np.pi/2, stepsizes)


# fpath = 'dipole_errors_0.01.npy'

"""
Below is how you save and load numpy arrays. If you're reading this, then you don't have
to run the dipole_error function again! The file will have pulled into your local respository
so you can just run this file and it will access the data straight away. :)
"""
### RUN THIS BLOCK TO GENERATE DATA AND SAVE IT ###
### DO THIS ONCE THEN ACCESS SAVED FILE TO SAVE TIME ###
# th_values, th_returns, deltas, lengths = dipole_error(50, -np.pi/2, np.pi/2, 0.01, 100000)
# with open(fpath, 'wb') as f:
#     np.save(f, [th_values, th_returns, deltas, lengths])

### RUN THIS BLOCK TO RETRIEVE SAVED DATA ###
# with open(fpath, 'rb') as f:
#     th_deltas = np.load(f, allow_pickle=True)
#     th_values, th_returns, deltas, lengths = th_deltas


#################### MEAN VALUES #######################

# th_values = [th for th in th_values if not np.isnan(th)]
# th_returns = [th for th in th_returns if not np.isnan(th)]
# deltas = [d for d in deltas if not np.isnan(d)]
# lengths = [l for l in lengths if not np.isnan(l)]
# th_gap = th_values[1]-th_values[0]
# print("theta gap=", th_gap)
# mean = np.mean(deltas)
# print("Mean Error (radians) = ", mean)
# mean_gap = mean/th_gap
# print(r"Mean Error / $\Delta\theta$ =", mean_gap)


#################### PLOTTING #######################

params = {
 'axes.labelsize': 14,
   'font.size': 14,
   'legend.fontsize': 14,
   'xtick.labelsize': 12,
   'ytick.labelsize': 12,
   'figure.figsize': [8,6]
   }
plt.rcParams.update(params)

#l = int(len(th_values)/2)

# fig, ax = plt.subplots(2,1, sharex=True)
# ax[0].plot(th_values[l:], deltas[l:]/th_gap, label="Step Size = 0.01")
# ax[0].plot(th_values[l:], [mean/th_gap for _ in th_values[l:]], label="Mean")
# ax[0].set_ylabel(r"(Angular Discrepancy)/$\Delta\theta$", fontsize = 'medium', labelpad = 17)
# ax[0].legend()
# ax[1].plot(th_values[l:], lengths[l:], label = "Step Size = 0.01")
# ax[1].set_ylabel("Fieldline Length (no. points)", fontsize='medium')
# ax[1].set_xlabel(r"$\theta$ (rad)", fontsize = 'medium')
# plt.legend()
# plt.rcParams.update(params)
# plt.show()

# plt.plot(th_values[l:], deltas[l:], label="Step Size = 0.01")
# plt.plot(th_values[l:], [mean for _ in th_values[l:]], label="Mean")
# plt.ylabel("Angular Discrepancy", fontsize = 'medium')
# plt.xlabel(r"$\theta$ (rad)", fontsize = 'medium')
# plt.legend()
# plt.show()

"""Dipole scaled errors vs latitude, with mean """
# plt.plot(th_values[l:], deltas[l:]/th_gap, label="Step Size = 0.01")
# plt.plot(th_values[l:], [mean_gap for _ in th_values[l:]], label="Mean")
# plt.ylabel(r"(Angular Discrepancy)/$\Delta\theta$", fontsize = 'medium')
# plt.xlabel(r"$\theta$ (rad)", fontsize = 'medium')
# plt.legend()
# plt.show()

# plt.rcParams.update(params)
# plt.plot(lengths[l:], deltas[l:]/th_gap, label = "Step Size = 0.01")
# plt.xscale('log')
# plt.xlabel("Log(Fieldline Length) [num. points]")
# plt.ylabel(r"(Angular Error)/$\Delta\theta$")
# plt.legend()
# plt.show()

### Want to observe uncertainty - use two field lines reflected in theta, look at difference between origin and ends.
# th_0 = 0.3
# x_1, y_1 = field_trace([1., th_0, 0.], dipole, 0.01, 100000)
# x_2, y_2 = field_trace([1., np.pi - th_0, 0.], dipole, 0.01, 100000)
# plt.plot(x_1, y_1, label="Northern Fieldline")
# plt.plot(x_2, y_2, label="Southern Fieldline")
# params = {
#  'axes.labelsize': 14,
#    'font.size': 18,
#    'legend.fontsize': 16,
#    'xtick.labelsize': 14,
#    'ytick.labelsize': 14,
#    'figure.figsize': [15,10]
#    }
# plt.rcParams.update(params)
# plt.xlabel("Distance in x (Scaled by Planetary Radius)")
# plt.ylabel("Distance in y (Scaled by Planetary Radius)")
# plt.legend()
# plt.show()


#################### ANALYTICAL VS PLOTTED #######################
def analytic_dipole_plot(numlines):
    theta_start = np.linspace(0, np.pi/2, numlines)
    def y(th, th_i):
        return np.sin(th)**2 * np.cos(th) / np.sin(th_i)**2
    def x(th, th_i):
        return np.sin(th)**3 / np.sin(th_i)**2

    for th_i in theta_start:
        coords = [(x(th, th_i), y(th, th_i)) for th in np.linspace(th_i, np.pi - th_i, 200)]
        x_arr, y_arr = map(list, zip(*coords))
        plt.plot(x_arr, y_arr, '-.', color = 'k')


# multiline_plot(25, th_max = np.pi/2)
# ls=(0, (3, 10, 1, 10, 1, 10))

        




