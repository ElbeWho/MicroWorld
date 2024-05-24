
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors 
import matplotlib.cbook as cbook
from matplotlib import cm
from mpl_toolkits import axes_grid1
from scipy.integrate import ode
from scipy.interpolate import RegularGridInterpolator as RGI

a=-4
b=4
steps = 0.01

mX, mY = np.mgrid[a:b:steps, a:b:steps]


u = 0.*mX
v = 0.*mY 
F = np.array([0, 1])
r0 = np.array([0, 0.001])
r=np.array([mX-r0[0], mY-r0[1]])
Id=np.array([[1,0],[0,1]])
Idf=np.dot(Id,F) 
rTF=(r*F[:,np.newaxis,np.newaxis]).sum(axis=0)
rrTF=(r*rTF[np.newaxis,])
modr=(r[0]**2+r[1]**2)**.5
u0, v0 =Idf[:,np.newaxis,np.newaxis]/modr+rrTF/modr**3.
u += u0
v += v0

print(mX[:, 0])
ui = RGI((mX[:, 0], mY[0, :]), u, method='linear')
vi = RGI((mX[:, 0], mY[0, :]), v, method='linear')

def Intergrate_Line(t, coord):
    # Function calculating Efield in every point for integrator
        global abort
        xi = coord[0]
        yi = coord[1]
        try:
            ex = ui([xi, yi])[0]
            ey = vi([xi, yi])[0]
            print(ey)
        except:
            abort = True
            return [False, False]
        n = (ex**2+ey**2)**0.5
        return [ex/n, ey/n]

integrmodel = 'vode' # 'lsoda', 'dopri5', 'dop853'
max_step = 2000
step_size = 0.01
fig, ax = plt.subplots()
dt = step_size


xstart = np.linspace(-3.9, 3.9, 21)

ystart = np.linspace(-3.9, -3.9, 21)
places=np.vstack([xstart,ystart]).T
        
print(places)

for p in places:

    r = ode(Intergrate_Line)
    r.set_integrator(integrmodel)
    lx=[p[0]]
    ly=[p[1]]
    r.set_initial_value([lx[0], ly[0]], 0)
    step = 0
    while r.successful():
                
        step += 1
        r.integrate(r.t+dt)
        x, y = r.y[0], r.y[1]

        lx.append(x)
        ly.append(y)
                
        if a >= x or b <= x or a >= y or b <= y or \
            step >= max_step:
            break
        if step >= max_step:
                    
            plt.plot(lx, ly, 'r') 
        else:
            plt.plot(lx, ly, 'k')
print("click me now")
plt.show()
"""if function == 'plot':
            self.steps = steps
            xx = np.linspace(self.a, self.b, steps)
            yy = np.linspace(self.a, self.b, steps)
            self.mX, self.mY = np.meshgrid(xx,yy)
            self.u = np.zeros(self.mX.shape)
            self.v = np.zeros(self.mY.shape)"""

xstart = np.linspace(-3.9, 3.9, 21)
ystart = np.linspace(-3.9, -3.9, 21)
r2 = np.array([2.0001, 2.0001])
f2 = np.array([0, -1])