import numpy as np
import Stokes_2D
r01 = np.array([0, 0.0001]) #defining positions of the force

F1 = np.array([0, -1]) #defining forces values
d = np.array([0, 1])
singularity1 = Stokes_2D.Equations2D(-4, 4, 0.01, "free") #creating an object
singularity1.dipole(r01, F1, d) #calling first stokeslet

singularity1.plot(direction=True, title="Dipole_puller.png") #plotting and saving
singularity1.show() #showing plot