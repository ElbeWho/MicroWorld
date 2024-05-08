import numpy as np
import Stokes_2D #importing software
import Stokes_3D 
#using class Stokes form MicroWorls
r0 = np.array([0, 0.001])
f = np.array([0, 1])
d = np.array([0, 0.5])
R = np.array([0,0,1])


plot_review = Stokes_2D.Equations2D(-4, 4, 0.08)
plot_review.stokeslet(r0, f)

plot_review.get_vi()
plot_review.streamlines()
plot_review.__show__()
