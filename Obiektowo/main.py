import numpy as np
import Stokes_2D 

r = np.array([0, 0.0001]) #defining position of the force
F = np.array([0, 1]) #defining force values
d = np.array([0, 1])

singularity1 = Stokes_2D.Equations2D(-4, 4, 0.01, "free") #creating an object

singularity1.dipole(r, d,-F) #calling stokeslet
singularity1.plot(direction=True, title="Dipole_pul.png") #plotting and saving
singularity1.show() #showing plot