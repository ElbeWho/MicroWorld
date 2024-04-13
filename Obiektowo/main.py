import numpy as np
import Stokes_2D #importing software
import Stokes_3D 
#using class Stokes form MicroWorls

r0 = np.array([0,0,0])
f = np.array([0, 1, 0])



plot = Stokes_3D.Equations3D(5, 5, 5) 
plot.int_params(1000, 0.01, 10, 0.1)
plot. stokeslet(r0, f)

plot.start()
plot.simulate()