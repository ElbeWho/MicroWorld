import numpy as np
import Stokes_2D 

r0 = np.array([0, 0.0001]) #defining position of the force
F = np.array([0, 1]) #defining force values

singularity1 = Stokes_2D.Equations2D(-4, 4, 0.01, "free") #creating an object
xstart = np.linspace(-3.9, 3.9, 21)
ystart = np.linspace(-3.9, -3.9, 21)
singularity1.stokeslet(r0, F) #calling stokeslet
singularity1.streamlines(xstart, ystart, export=True) #plotting and saving
singularity1.show() #showing plot
