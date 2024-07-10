import numpy as np
import Stokes_2D
r01 = np.array([0, 0.5001]) #defining positions of the force
F1 = np.array([1, 1]) #defining forces values
r02 = np.array([3.0001, 3.0001]) #defining positions of the force
F2 = np.array([2, 2]) #defining forces values
r03 = np.array([-3, 2.0001]) #defining positions of the force
F3 = np.array([0.5, 3])
singularity1 = Stokes_2D.Equations2D(-4, 4, 0.01, "free") #creating an object
singularity1.totality_per(r01, F1, 3) #ok

singularity1.plot() #plotting and saving
singularity1.show() #showing plot