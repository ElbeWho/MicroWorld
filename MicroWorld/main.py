import numpy as np
import Stokes_2D
r01 = np.array([-2, 2.0001]) #defining positions of the force
F = np.array([1, 1])
r02 = np.array([1.5, 2.0001])
r03 = np.array([3, 2.0001])
singularity1 = Stokes_2D.Equations2D(-4, 4, 0.01, "wall") #creating an object
singularity1.totality_par(r01, -F, 3)
singularity1.totality_par(r02, F, 3)
singularity1.totality_par(r03, -F, 3)
singularity1.plot(title="example2_par_wall.png")
singularity1.show() #showing plot