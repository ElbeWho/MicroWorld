import numpy as np
import Stokes_2D 

r0 = np.array([0, .5]) #defining positions of the force
d = np.array([0, 1])
F = np.array([1, 1]) #defining forces values
Frel = np.array([0, 1])

singularity1 = Stokes_2D.Equations2D(-4, 4, 0.01, "free") #creating an object

singularity1.hard_wall_par(r0, F) #calling first stokeslet


singularity1.plot(title="dip.png") #plotting and saving
singularity1.show()
