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

#monopol.source_doublet(r0, f)
#nie no dodwanie tego jest kompletnie bez sensu -.-
#najlepiej, eby prędkość się aktualizowała przy kadym wywołaniu funkcji 
#troche głupie to ale wsm ma sens bo takie jest przeznaczenie tego modułu




monopol.stokeslet(f, r0)


monopol1 = custom_functions.Stokes("Monopol1") 

monopol1.__meshgrid__(5, 100)

monopol1.stokeslet(-f, r1)

dipol = monopol + monopol1

dipol.__plot__()
dipol.__show__()