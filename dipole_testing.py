"""
PyFields dipole_testing.py

28.10.2020

Testing model using a dipole
"""

from test_wrapper import *
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pandas as pd

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

th_values = np.linspace((-np.pi)/2, (np.pi)/2, 50)
for th in th_values:
    if th==0 or th==np.pi or th==2*np.pi:
        pass
    else:
         x, y = field_trace([1., th, 0.1], dipole, 0.01, 100000)
         plt.plot(x, y, color = 'b')
         print(th)

plt.show()
