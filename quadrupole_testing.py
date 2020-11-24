from all_funcs import *
import matplotlib.pyplot as plt
import matplotlib as mat
import numpy as np
from collections import Counter


########################## QUADRUPOLE TEST ###############################

params = {
 'axes.labelsize': 14,
   'font.size': 14,
   'legend.fontsize': 14,
   'xtick.labelsize': 12,
   'ytick.labelsize': 12,
   'figure.figsize': [8,6]
   }
plt.rcParams.update(params)
multiline_plot(50, coeffs=quadrupole)
plt.xlabel("Distance in x")
plt.ylabel("Distance in y")
#plt.annotate("Traced Quadropole, ds=0.01", xy=(1,1))
plt.legend((mat.lines.Line2D([0,0], [1,1], color = 'r'),mat.lines.Line2D([0,0], [1,1], color = 'b')), ('Southbound Line','Northbound Line'))
plt.show()


########################## ANALYTICAL FIELD COMPARISON ##########################

params = {
 'axes.labelsize': 14,
   'font.size': 14,
   'legend.fontsize': 14,
   'xtick.labelsize': 12,
   'ytick.labelsize': 12,
   'figure.figsize': [8,6]
   }
plt.rcParams.update(params)

field = 'Quadrupole'

multiline_plot(25, th_max = np.pi/2, coeffs = quadrupole)
_analytic_field_plot(0, np.pi/2, 25, 0.001, 'quad')
plt.legend((mat.lines.Line2D([0,0], [1,1], color = 'r'),mat.lines.Line2D([0,0], [1,1], color = 'k', linestyle='--')), (f'Traced {field}, ds = 0.01',f'Analytical {field}'))
plt.show()
