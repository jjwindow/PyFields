"""
PyFields dipole_testing.py

28.10.2020

Testing model using a dipole
"""

from all_funcs import field_trace, dipole
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm

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
# plt.show()

# Want to observe uncertainty - use two field lines reflected in theta, look at difference between origin and ends.
# th_0 = 0.3
# x_1, y_1 = field_trace([1., th_0, 0.], dipole, 0.01, 100000)
# x_2, y_2 = field_trace([1., np.pi - th_0, 0.], dipole, 0.01, 100000)
# plt.plot(x_1, y_1, label="Northern fieldline")
# plt.plot(x_2, y_2, label="Southern fieldline")
# plt.xlabel("x scaled by planetary radius")
# plt.ylabel("y scaled by planetary radius")
# plt.legend()
# plt.show()