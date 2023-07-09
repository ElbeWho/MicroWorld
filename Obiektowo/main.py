import numpy as np
import ast
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors 
import matplotlib.cbook as cbook

from matplotlib import cm
import sympy

from mpl_toolkits import axes_grid1

import custom_functions as cf

pierwszy = cf.Stokes(5, 100)

r0=np.array([0,0])            # position of the force red arrow - touchdown point
f=np.array([0,1])             # direction of the force

pierwszy.stokeslet(f)

