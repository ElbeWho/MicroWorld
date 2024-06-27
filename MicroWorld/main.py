import numpy as np
import Stokes_2D 

r01 = np.array([0, -.50001]) #defining positions of the force
r02 = np.array([-1, .50001])
r03 = np.array([1, .50001])
F1 = np.array([0, 1]) #defining forces values
F2 = np.array([0, -.5])
F3 = np.array([0, -.5])
singularity1 = Stokes_2D.Equations2D(-4, 4, 0.01, "free") #creating an object

singularity1.stokeslet(r01, F1) #calling first stokeslet
singularity1.stokeslet(r02, F2) #calling second stokeslet
singularity1.stokeslet(r03, F3) #calling second stokeslet
singularity1.plot(direction=True, title="alga_mean.png") #plotting and saving
singularity1.show() #showing plot