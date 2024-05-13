import numpy as np
import Stokes_2D #importing software
import Stokes_3D 
#using class Stokes form MicroWorls
r0 = np.array([0, 1.001])
f = np.array([0, 1])
d = np.array([0, 0.5])
R = np.array([0,0,1])


plot_review = Stokes_2D.Equations2D(-5, 5, 0.01, "wall")

plot_review.stokes_dipole(r0, f)
#plot_review.hard_wall_par(r0, f)



plot_review.__plot__()
plot_review.__show__()
#plot_review.save_plot("wall_first2")
