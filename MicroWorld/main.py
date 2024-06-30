import numpy as np
import Stokes_2D
r01 = np.array([0, 0.5001]) #defining positions of the force

F1 = np.array([1, 1]) #defining forces values
d = np.array([0, 1])
singularity1 = Stokes_2D.Equations2D(-4, 4, 0.01, "free") #creating an object
singularity1.totality_par(r01, F1, 1) #calling first stokeslet

singularity1.plot( title="jeedynka_par.png") #plotting and saving
singularity1.show() #showing plot