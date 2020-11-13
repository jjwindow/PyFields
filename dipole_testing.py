"""
PyFields dipole_testing.py

28.10.2020

Testing model using a dipole
"""

from all_funcs import *
import matplotlib.pyplot as plt
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
    for i, th in enumerate(th_values):
        start_pos = [1., th, 0.]
        field = field_trace(start_pos, dipole, ds, max_iter, axes=None)
        if field is not None:
            (p_arr, B_arr) = field
            th_final = p_arr[-1][1]
            deltas.append(abs(th_final-th_returns[i]))
        else:
            th_values[i] = np.nan
    th_values = [th for th in th_values if not np.isnan(th)]    
    deltas = np.array(deltas)
      
    # print("angle at which dipole lines leave surface=", th_values)
    # print("angle at which dipole lines enter surface=", th_returns)
    # print("angular discrepancies=", deltas)
    return th_values, th_returns, deltas

fpath = 'dipole_errors_0.01.npy'

"""
Below is how you save and load numpy arrays. If you're reading this, then you don't have
to run the dipole_error function again! The file will have pulled into your local respository
so you can just run this file and it will access the data straight away. :)
"""
### RUN THIS BLOCK TO GENERATE DATA AND SAVE IT ###
### DO THIS ONCE THEN ACCESS SAVED FILE TO SAVE TIME ###
# th_values, th_returns, deltas = dipole_error(50, -np.pi/2, np.pi/2, 0.01, 100000)
# with open(fpath, 'wb') as f:
#     np.save(f, [th_values, deltas])

### RUN THIS BLOCK TO RETRIEVE SAVED DATA ###
with open(fpath, 'rb') as f:
    th_deltas = np.load(f, allow_pickle=True)
    th_values, deltas = th_deltas[0], th_deltas[1]

th_gap = th_values[1]-th_values[0]
mean = deltas.mean()
print("Mean Error (radians) = ", mean)

plt.plot(th_values, deltas/th_gap, label="Step Size = 0.01")
plt.plot(th_values, [mean/th_gap for _ in th_values], label="Mean")
params = {
 'axes.labelsize': 16,
   'font.size': 14,
   'legend.fontsize': 14,
   'xtick.labelsize': 16,
   'ytick.labelsize': 16,
   'figure.figsize': [15,10]
   }
plt.rcParams.update(params)
plt.xlabel(r"$\theta$ (rad)", fontsize = 'medium')
plt.ylabel(r"(Angular Discrepancy)/$\Delta\theta$", fontsize = 'medium')
plt.legend()
plt.show()


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





