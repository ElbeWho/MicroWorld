import numpy as np
import MicroWorld #importing software
#using class Stokes form MicroWorls
plot_for_posem = MicroWorld.Stokes(5, 100) 
#radiuses and distance d
r0 = np.array([0,0])
r1 = np.array([2, 2])
r2 = np.array([-2, -2])
d = np.array([0,0.1])
#forces        
f  = np.array([0.55, 0.3]) 
#adding choosen characteristic sollutions
plot_for_posem.source(r0)
plot_for_posem.dipole(r1, d, f)
plot_for_posem.stokeslet(f, r2)
plot_for_posem.__plot__()
plot_for_posem.__show__()

