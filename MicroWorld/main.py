import numpy as np
import Stokes_2D 

r01 = np.array([0, 0.0001]) #defining positions of the force
r02 = np.array([0, -0.1001])
F = np.array([0, 1])
F1 = np.array([1, 0]) #defining forces values
F2 = np.array([0, -3])

singularity1 = Stokes_2D.Equations2D(-4, 4, 0.01, "free") #creating an object
P=1
singularity1.source(r01, P) #calling first stokeslet
singularity1.plot( title="source.png") #plotting and saving
singularity1.show() #showing plot