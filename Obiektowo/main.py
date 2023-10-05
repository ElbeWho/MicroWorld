import numpy as np
import ast
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors 
import matplotlib.cbook as cbook
from matplotlib import cm
import sympy
from mpl_toolkits import axes_grid1

import custom_functions

monopol = custom_functions.Stokes("Monopol")
monopol.__meshgrid__(5, 100)
r0=np.array([0.0,0.0])            # position of the force red arrow - touchdown point
f0=np.array([0,1])             # direction of the force
r1=np.array([0.0,0.0])            # position of the force red arrow - touchdown point
f1=np.array([0,-1])   
r2=np.array([1,0])            # position of the force red arrow - touchdown point
f2=np.array([0,-1]) 
r3=np.array([4,0])            # position of the force red arrow - touchdown point
f3=np.array([0,2])


monopol.forces(2)          # direction of the force
monopol.entries([ [r1,f1],  [r2,f2], [r3, f3]])
monopol.many_stokeslets()
monopol.get_entries()
#monopol.stokeslet(f0, r0)

monopol.__plot__()
monopol.save_plot()