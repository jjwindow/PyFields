"""
PyFields dipole_testing.py

28.10.2020

Testing model using a dipole
"""

from all_funcs import *
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

### PRINTING INDIVIDUAL LINES ###

#pstart_1 = [1., 1.67, 0.1]
# pstart_2 = [1., 0.15, 0.1]
# pstart_3 = [1., 0.2, 0.1]
# pstart_4 = [1., 0.25, 0.1]

#x_1, y_1 = field_trace(pstart_1, dipole, 0.01, 100000)
# x_2, y_2 = field_trace(pstart_2, dipole, 0.01, 100000)
# x_3, y_3 = field_trace(pstart_3, dipole, 0.01, 100000)
# x_4, y_4 = field_trace(pstart_4, dipole, 0.01, 100000)


#plt.plot(x_1, y_1, color = 'b')
# plt.plot(x_2, y_2, color = 'b')
# plt.plot(x_3, y_3, color = 'b')
# plt.plot(x_4, y_4, color = 'b')
# plt.title("Test Dipole")
#plt.show()


# multilines(50)
# plt.title("Test Dipole")
# plt.show()

### Want to observe uncertainty - use two field lines reflected in theta, look at difference between origin and ends.
# th_0 = 0.3
# x_1, y_1 = field_trace([1., th_0, 0.], dipole, 0.01, 100000)
# x_2, y_2 = field_trace([1., np.pi - th_0, 0.], dipole, 0.01, 100000)
# plt.plot(x_1, y_1, label="Northern fieldline")
# plt.plot(x_2, y_2, label="Southern fieldline")
# plt.xlabel("x scaled by planetary radius")
# plt.ylabel("y scaled by planetary radius")
# plt.legend()
# plt.show()


####### ERROR CALCULATION ######

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

th_values, th_returns, deltas = dipole_error(50, -np.pi/2, np.pi/2, 0.01, 100000)
th_gap = th_values[1]-th_values[0]
mean = deltas.mean()

plt.plot(th_values, deltas/th_gap, label="Step Size = 0.01")
plt.plot(th_values, [mean/th_gap for _ in th_values], label="Mean")
plt.xlabel(r"$\theta$ (rad)")
plt.ylabel(r"(Angular Discrepancy)/$\Delta\theta$")
plt.legend()
plt.show()



# th_values = np.linspace(-np.pi/2, np.pi/2, 50)
# th_returns = []
# for th in th_values:
#     if (th < 0.0):
#         th_return = -np.pi + abs(th)
#         th_returns.append(th_return)
#     if (th > 0.0):
#         th_return = np.pi - th
#         th_returns.append(th_return)
# th_returns = np.array(th_returns)
# for th in th_values:
#     start_pos = [1., th, 0.]
#     field_trace(start_pos, dipole, ds, max_iter, axes=None)
    
# print("angle at which dipole lines leave surface=", th_values)
# print("angle at which dipole lines enter surface=", th_returns)

