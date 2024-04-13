
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import ode
from scipy.interpolate import RegularGridInterpolator as RGI

class Equations3D:

    def __init__ (self, a: float, b: float, c: float):

        self.a = a
        self.b = b 
        self.c = c

        self.X, self.Y, self.Z = np.mgrid[a:b:c, a:b:c, a:b:c]

        self.U = 0.*self.X
        self.V= 0.*self.Y
        self.W = 0.*self.Z


    def int_params(self, max_step, step_size, n, Rsphere):

        self.integrmodel = 'vode' # 'lsoda', 'dopri5', 'dop853'
        self.max_step = max_step
        self.step_size = step_size
        self.n_sphere = n
        self.Rsphere = Rsphere
        self.R2sphere = Rsphere**2

    def stokeslet(self, r0, f):

        Id = np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        r = np.array([self.X - r0[0], self.Y - r0[1], self.Z - r0[2]])

        Idf = np.dot(Id, f)
        rTf = (r * f[:, np.newaxis, np.newaxis, np.newaxis]).sum(axis=0)
        rrTf = (r * rTf[np.newaxis,])
        modr = np.sqrt(r[0]**2 + r[1]**2 + r[2]**2)

        u, v, w = Idf[:, np.newaxis, np.newaxis, np.newaxis] / modr + rrTf / modr**3.

        self.U += u
        self.V += v
        self.W += w

        
        self.Ui = RGI((self.X[:, 0, 0], self.Y[0, :, 0], self.Z[0, 0, :]), u, method='linear')
        self.Vi = RGI((self.X[:, 0, 0], self.Y[0, :, 0], self.Z[0, 0, :]), v, method='linear')
        self.Wi = RGI((self.X[:, 0, 0], self.Y[0, :, 0], self.Z[0, 0, :]), w, method='linear') 
    
    def IL3(self, t, coord):
        xi, yi, zi = coord
        try:
            ex = self.Ui([xi, yi, zi])[0]
            ey = self.Vi([xi, yi, zi])[0]
            ez = self.Wi([xi, yi, zi])[0]
        except:
            return [0, 0, 0]
        return [ex, ey, ez]


    def start(self):
        self.n = 5

        self.space = np.linspace(-2.5, 2.5, self.n)

        self.exes = np.zeros(self.n*self.n)
        self.ezes = np.zeros(self.n*self.n)
        self.eyes = np.zeros(self.n*self.n)

        self.step = 0

        self.exes = []
        self.ezes = []
        self.eyes = []

        for i in range(0, self.n):
            for k in range(0, self.n):
                self.step = self.step +1

            self.exes.append(self.space[i])
            self.ezes.append(self.space[k])
            self.eyes.append(-2.5)

        self.xxx = np.zeros(self.n*self.n)
        self.yyy = np.zeros(self.n*self.n)
        self.zzz = np.zeros(self.n*self.n)

        for i in range(0, self.n*self.n):
            self.xxx[i] = self.exes[i]
            self.yyy[i] = self.eyes[i]
            self.zzz[i] = self.ezes[i]
        
        self.places = np.vstack([self.xxx, self.yyy, self.zzz]).T

    def simulate(self):

        self.ax = plt.figure().add_subplot(projection='3d')

        dt = self.step_size

        self.xlist = []
        self.ylist = []
        self.zlist = []
        for p in self.places:

            r = ode(self.IL3)
            r.set_integrator(self.integrmodel)
            lx=[p[0]]
            ly=[p[1]]
            lz = [p[2]]
            r.set_initial_value([lx[0], ly[0], lz[0]], 0)

            step = 0
            abort = False

        while r.successful():
            self.step += 1

            r.integrate(r.t+dt)
            x, y, z = r.y[0], r.y[1], r.y[2]

            lx.append(x)
            ly.append(y)
            lz.append(z)
            self.xlist.append(x)
            self.ylist.append(y)
            self.zlist.append(z)
      
            if self.s >= x or self.e <= x or self.s >= y or self.e <= y or  self.s >= z or self.e <= z:
                self.xlist.append(None)
                self.ylist.append(None)
                self.zlist.append(None)
                break

            if step >= self.max_step:
                self.xlist.append(None)
                self.ylist.append(None)
                self.zlist.append(None) 
                break

            if abort is True:
                abort = False
                self.xlist.append(None)
                self.ylist.append(None)
                self.zlist.append(None) 
                break

        plt.plot(lx, ly, lz,  'k', lw=0.5)

        plt.show()
