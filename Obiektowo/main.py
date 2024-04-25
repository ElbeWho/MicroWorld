import numpy as np
import Stokes_2D #importing software
import Stokes_3D 
#using class Stokes form MicroWorls
r0 = np.array([0,0])
f = np.array([0, 1])
d = np.array([0, 0.5])
R = np.array([0,0,1])
plot = Stokes_3D.Equations3D(-4, 4, 0.05) 


plot_review = Stokes_2D.Equations2D(4, 4, 100)
plot_review.rotlet_R(r0, R)
plot_review.source_dipole()

plot_review.__plot__()
plot_review.__show__()
plot_review.save_plot('pierwszy')