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
########################## DIPOLE TEST ###############################

# multilines(50)
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
# plt.show()


########################## ERROR CALCULATION ##########################

def dipole_error(num, th_min, th_max, ds, max_iter):
    th_values = np.linspace(th_min, th_max, num)
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
    with tqdm(total=len(th_values), desc="Tracing Fields...") as bar:
        for i, th in enumerate(th_values):
            start_pos = [1., th, 0.]
            field = field_trace(start_pos, dipole, ds, max_iter, axes=None)
            if field is not None:
                (p_arr, B_arr) = field
                th_final = p_arr[-1][1]
                lengths.append(len(p_arr))
                deltas.append(abs(th_final-th_returns[i]))
            else:
                th_values[i] = np.nan
            bar.update()
    th_values = [th for th in th_values if not np.isnan(th)]    
    deltas = np.array(deltas)
    lengths = np.array(lengths)

    return th_values, th_returns, deltas, lengths

fpath = 'dipole_errors_0.01.npy'

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
with open(fpath, 'rb') as f:
    th_deltas = np.load(f, allow_pickle=True)
    th_values, th_returns, deltas, lengths = th_deltas

th_gap = th_values[1]-th_values[0]
mean = deltas.mean()
print("Mean Error (radians) = ", mean)

#################### PLOTTING #######################

params = {
 'axes.labelsize': 14,
   'font.size': 14,
   'legend.fontsize': 14,
   'xtick.labelsize': 12,
   'ytick.labelsize': 12,
   'figure.figsize': [8,6]
   }
l = int(len(th_values)/2)

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

plt.plot(th_values[l:], deltas[l:]/th_gap, label="Step Size = 0.01")
plt.plot(th_values[l:], [mean/th_gap for _ in th_values[l:]], label="Mean")
plt.ylabel(r"(Angular Discrepancy)/$\Delta\theta$", fontsize = 'medium')
plt.xlabel(r"$\theta$ (rad)", fontsize = 'medium')
plt.legend()
plt.show()

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

def analytic_dipole_plot(numlines):
    theta_start = np.linspace(0, np.pi/2, numlines)
    def y(th, th_i):
        return np.sin(th)**2 * np.cos(th) / np.sin(th_i)**2
    def x(th, th_i):
        return np.sin(th)**3 / np.sin(th_i)**2

    for th_i in theta_start[1:]:
        coords = [(x(th, th_i), y(th, th_i)) for th in np.linspace(th_i, np.pi - th_i, 200)]
        x_arr, y_arr = map(list, zip(*coords))
        plt.plot(x_arr, y_arr, '-.', color = 'purple')
    

# multiline_plot(25, th_max = np.pi/2)
# analytic_dipole_plot(25)
# plt.legend((mat.lines.Line2D([0,0], [1,1], color = 'r'),mat.lines.Line2D([0,0], [1,1], color = 'purple')), ('Traced Dipole, ds = 0.01','Analytical Dipole'))
# plt.show()
        




