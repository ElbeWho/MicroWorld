import numpy as np
import custom_functions

monopol = custom_functions.Stokes(5, 100) 
                                  
'''r'$\displaystyle\\v(r)='
               r'\frac{\mathbf{F}}{8 \pi \eta r } ( \mathds{1} + \frac{\mathbf{rr}}{r^2})$' '''


r0 = np.array([0, 0.5])
r2 = np.array([0, -0.5])
r1 = np.array([0, 3.5])  
f1 = np.array([1, 3.5])          

f = np.array([0,1])   

#monopol.source_doublet(r0, f)
#nie no dodwanie tego jest kompletnie bez sensu -.-
#najlepiej, zeby prędkość się aktualizowała przy kadym wywołaniu funkcji 
#troche głupie to ale wsm ma sens bo takie jest przeznaczenie tego modułu

#lepiej skala nie na sztywno
#ale przydałoby móc się zablokować na jakichś wartościach


monopol.stokeslet(f, r0)

monopol.source_doublet(-f, r1)



monopol.__plot__()
monopol.__show__()