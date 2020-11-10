from all_funcs import field_trace
import matplotlib.pyplot as plt
import numpy as np

# Quadrupole coefficients
g_Q = np.array([[0., 0., 0., 0.], [0., 0., 0., 0.], [1., 0., 0., 0.], [0., 0., 0., 0.]])
h_Q = np.array([[0., 0., 0., 0.], [0., 0., 0., 0.], [1., 0., 0., 0.], [0., 0., 0., 0.]])
a_Q = 1

quadrupole = (a_Q, g_Q, h_Q)

# def multilines(num):
#      th_values = np.linspace(0, 2*(np.pi), num)
#      for th in th_values:
#          if th==0 or th==np.pi or th==2*np.pi:
#              pass
#          else:
#              x, y = field_trace([1., th, 0.1], quadrupole, 0.01, 100000)
#              plt.plot(x, y, color = 'b')

# multilines(100)
# plt.title("Test Quadrupole")
# plt.show()

th_values = np.linspace(0, 2*(np.pi), 50)
print(th_values)

a=np.empty(50)
a.fill(np.pi)
print(a)

b = a - th_values
print(b)