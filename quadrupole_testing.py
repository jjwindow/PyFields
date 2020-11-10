from all_funcs import field_trace
from dipole_testing import multilines
import matplotlib.pyplot as plt
import matplotlib as mat
import numpy as np
from collections import Counter


# Quadrupole coefficients
g_Q = np.array([[0., 0., 0., 0.], [0., 0., 0., 0.], [1., 0., 0., 0.], [0., 0., 0., 0.]])
h_Q = np.array([[0., 0., 0., 0.], [0., 0., 0., 0.], [1., 0., 0., 0.], [0., 0., 0., 0.]])
a_Q = 1

quadrupole = (a_Q, g_Q, h_Q)

multilines(25, th_min = 0, th_max=np.pi/4, coeffs=quadrupole)
multilines(25, th_min = 3*np.pi/4, th_max=np.pi, coeffs=quadrupole)
multilines(25, th_min = np.pi, th_max=5*np.pi/4, coeffs=quadrupole)
multilines(25, th_min=7*np.pi/4, th_max = 2*np.pi, coeffs = quadrupole)
plt.title("Test Quadrupole")
plt.legend((mat.lines.Line2D([0,0], [1,1], color = 'r'),mat.lines.Line2D([0,0], [1,1], color = 'b')), ('Southbound Line','Northbound Line'))
plt.show()


#### THETA SYMMETRY TESTING #####
# th_values = np.linspace(0, 2*(np.pi), 50)
# print(th_values)

# a=np.empty(50)
# a.fill(np.pi)
# print(a)

# b = abs(a - th_values)
# print(b)
# precision = 6
# print(Counter(map(round, b, [precision for _ in b])))