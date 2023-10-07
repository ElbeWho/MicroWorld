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

r0 = np.array([0, 0])
r1 = np.array([0, 3.5])            
f = np.array([0,1])   
f1 = np.array([1, 3.5])
       
monopol.entries([ [f,r0], [f, r1]])

monopol.many_stokeslets()

monopol.__plot__()

#monopol.save_plot()