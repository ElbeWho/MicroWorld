import numpy as np
import MicroWorld #importing software
#using class Stokes form MicroWorls
plot = MicroWorld.Stokes(5, 100) 
#radiuses and distance d
r0 = np.array([0,0])

#adding choosen characteristic solutions
R = np.array([1, 0])
plot.Stokes_dipole( R, r0)

plot.__plot__()
plot.__show__()

