import numpy as np
import Stokes_2D
r01 = np.array([0, -0.5001]) #defining positions of the force
F1 = np.array([0, 1]) #defining forces values
r02 = np.array([-1, 0.5001])  #defining positions of the force
F2 = np.array([0, -0.5]) #defining forces values
r03 = np.array([1, 0.5001])  #defining positions of the force
F3 = np.array([0, -0.5])
singularity1 = Stokes_2D.Equations2D(-10, 10, 0.02, "free") #creating an object
singularity1.stokeslet(r01, F1)#ok
singularity1.stokeslet(r02, F2)#ok
singularity1.stokeslet(r03, F3)#ok
singularity1.plot(direction=True, title="alga_mean_dipole.png") #plotting and saving
singularity1.show() #showing plot