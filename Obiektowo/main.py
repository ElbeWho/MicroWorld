import numpy as np
import custom_functions

monopol = custom_functions.Stokes("Monopol") 
                                  
'''r'$\displaystyle\\v(r)='
               r'\frac{\mathbf{F}}{8 \pi \eta r } ( \mathds{1} + \frac{\mathbf{rr}}{r^2})$' '''

monopol.__meshgrid__(5, 100)

r0 = np.array([0, 0.5])
r2 = np.array([0, -0.5])
r1 = np.array([0, 3.5])  
f1 = np.array([1, 3.5])          

f = np.array([0,1])   



d = np.array([1,0])

monopol.stresslet(f, r0, d)
 
#monopol.entries([ [f,r0], [f1, r0]])

#monopol.many_stokeslets()

monopol.__plot__()

monopol.__show__()

monopol.title()

#monopol.save_plot()