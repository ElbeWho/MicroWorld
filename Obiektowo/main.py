import numpy as np
import Stokes_2D #importing software
import Stokes_3D 
#using class Stokes form MicroWorls
r0 = np.array([0, 0.0001])
r2 = np.array([-1, 0])
f = np.array([0, 1])
d = np.array([0, 0.5])
R = np.array([0,0,1])


plot_review = Stokes_2D.Equations2D(-4, 4, 0.01, "free")


xstart = np.linspace(-3.9, 3.9, 21)
ystart = np.linspace(-3.9, -3.9, 21)
plot_review.stokeslet(r0, f)
plot_review.streamlines(xstart, ystart, mesh=True)
plot_review.__show__() 
