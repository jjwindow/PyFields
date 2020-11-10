"""
PyFields dipole_testing.py

28.10.2020

Testing model using a dipole
"""

from all_funcs import field_trace, dipole
import matplotlib.pyplot as plt
import numpy as np

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

def multilines(num, th_min = 0, th_max = np.pi, coeffs = dipole, ds = 0.01, maxits = 100000):
    """
    Plots 'num' (int) field lines for equally spaced theta values between th_min and th_max.
    Field lines calculated using field coefficients given by coeffs (tuple), stepsize ds (float),
    and terminating after maxits (int). Use plt.show() to display plot after calling.

    MOVE TO ALL_FUNCS EVENTUALLY
    """
    th_values = np.linspace(th_min, th_max, num)
    for th in th_values:
        if th==0 or th==np.pi or th==2*np.pi:
            pass
        else:
            x, y = field_trace([1., th, 0.], coeffs, ds, maxits)
            if y[0] > y[-1]:
                colour = 'r'
            else:
                colour = 'b'
            plt.plot(x, y, color = colour)
            #  print(th)

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