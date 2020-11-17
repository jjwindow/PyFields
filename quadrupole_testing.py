from all_funcs import *
import matplotlib.pyplot as plt
import matplotlib as mat
import numpy as np
from collections import Counter


# Quadrupole coefficients
g_Q = np.array([[0., 0., 0., 0.], [0., 0., 0., 0.], [1., 0., 0., 0.], [0., 0., 0., 0.]])
h_Q = np.array([[0., 0., 0., 0.], [0., 0., 0., 0.], [1., 0., 0., 0.], [0., 0., 0., 0.]])
a_Q = 1

quadrupole = (a_Q, g_Q, h_Q)

########################## QUADRUPOLE TEST ###############################

# params = {
#  'axes.labelsize': 14,
#    'font.size': 14,
#    'legend.fontsize': 14,
#    'xtick.labelsize': 12,
#    'ytick.labelsize': 12,
#    'figure.figsize': [8,6]
#    }
# plt.rcParams.update(params)
# # multilines(25, th_min = 0, th_max=np.pi/4, coeffs=quadrupole)
# # multilines(25, th_min = 3*np.pi/4, th_max=np.pi, coeffs=quadrupole)
# # multilines(25, th_min = np.pi, th_max=5*np.pi/4, coeffs=quadrupole)
# # multilines(25, th_min=7*np.pi/4, th_max = 2*np.pi, coeffs = quadrupole)
# multiline_plot(100, th_min=0.0, th_max=2*np.pi, coeffs=quadrupole)
# plt.xlabel("Distance in x (Scaled by Planetary Radius)")
# plt.ylabel("Distance in y (Scaled by Planetary Radius)")
# plt.legend((mat.lines.Line2D([0,0], [1,1], color = 'r'),mat.lines.Line2D([0,0], [1,1], color = 'b')), ('Southbound Line','Northbound Line'))
# plt.show()


########################## THETA SYMMETRY TESTING ##########################
# th_values = np.linspace(0, 2*(np.pi), 50)
# print(th_values)

# a=np.empty(50)
# a.fill(np.pi)

# b = abs(a - th_values)
# print(b)

# precision = 6
# print(Counter(map(round, b, [precision for _ in b])))


########################## ERROR CALCULATION ##########################

def quad_error(num, th_min, th_max, ds, max_iter):
    th_values = np.linspace(th_min, th_max, num)
    th_returns = []
    for th in th_values:
        if ((th >= 0.0) and  ( th <= np.pi/2)):
            th_return = np.pi/2 - th
            th_returns.append(th_return)
        elif ((th > np.pi/2) and (th <= np.pi)):
            th_return = (3*np.pi)/2 - th
            th_returns.append(th_return)
        elif ((th > np.pi) and (th <= (3*np.pi)/2)):
            th_return = (5*np.pi)/2 - th
            th_returns.append(th_return)
        else:
            th_return = (7*np.pi)/2 - th
            th_returns.append(th_return)
    
    th_returns = np.array(th_returns)

    deltas = []
    lengths = []

    
    for i, th in enumerate(th_values):
        if th==0 or th==np.pi/2 or th==np.pi or th==(3*np.pi)/2 or th==2*np.pi:
            deltas.append(np.nan)
            lengths.append(np.nan)
        else:
            start_pos = [1., th, 0.]
            field = field_trace(start_pos, quadrupole, ds, max_iter, axes=None)
            if field is not None:
                (p_arr, B_arr) = field
                th_final = p_arr[-1][1]
                lengths.append(len(p_arr))
                deltas.append(abs(th_final-th_returns[i]))
            else:
                th_values[i] = np.nan
    th_values = [th for th in th_values if not np.isnan(th)]    
    deltas = np.array([deltas)
    lengths = np.array(lengths)

    return th_values, th_returns, deltas, lengths

    fpath = 'quadrupole_errors_0.01.npy'

"""
Below is how you save and load numpy arrays. If you're reading this, then you don't have
to run the dipole_error function again! The file will have pulled into your local respository
so you can just run this file and it will access the data straight away. :)
"""
### RUN THIS BLOCK TO GENERATE DATA AND SAVE IT ###
### DO THIS ONCE THEN ACCESS SAVED FILE TO SAVE TIME ###
th_values, th_returns, deltas, lengths = quad_error(50, 0, 2*np.pi, 0.01, 100000)
with open(fpath, 'wb') as f:
    np.save(f, [th_values, th_returns, deltas, lengths])

### RUN THIS BLOCK TO RETRIEVE SAVED DATA ###
with open(fpath, 'rb') as f:
    th_deltas = np.load(f, allow_pickle=True)
    th_values, th_returns, deltas, lengths = th_deltas

th_gap = th_values[1]-th_values[0]
mean = deltas.mean()
print("Mean Error (radians) = ", mean)
mean_gap = mean/th_gap
print(r"Mean Error / $\Delta\theta$ =", mean_gap)


params = {
 'axes.labelsize': 14,
   'font.size': 14,
   'legend.fontsize': 14,
   'xtick.labelsize': 12,
   'ytick.labelsize': 12,
   'figure.figsize': [8,6]
   }
