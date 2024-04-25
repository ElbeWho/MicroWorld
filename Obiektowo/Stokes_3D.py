
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
        self.U = 0*self.X
        self.V= 0*self.Y
        self.W = 0*self.Z
        


    def int_params(self, max_step, step_size, z, Rsphere):

        self.integrmodel = 'vode' # 'lsoda', 'dopri5', 'dop853'
        self.max_step = max_step
        self.step_size = step_size
        self.n_sphere = z
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
        
    def rotlet(self, r0, d, e):

        r = np.array([self.X - r0[0], self.Y - r0[1], self.Z - r0[2]])
        modr = (r[0]**2+r[1]**2+r[2]**2)**.5 
        jeden = ((d[0]*r[0]+ d[1]*r[1]+d[2]*r[1])*e[:, np.newaxis, np.newaxis, np.newaxis] 
                 - (e[0]*r[0]+e[1]*r[1]+e[2]*e[2])*d[:, np.newaxis, np.newaxis, np.newaxis] )/modr**3
        u, v, w = jeden
        self.U += u
        self.V += v
        self.W += w

    
    
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

        self.exes = []
        self.ezes = []
        self.eyes = []

        self.step = 0

        for i in range(0, self.n):
            for k in range(0, self.n):
                self.step = self.step + 1
                self.exes.append(self.space[i])
                self.ezes.append(self.space[k])
                self.eyes.append(-2.5)

        self.xxx = np.zeros(self.n * self.n)
        self.yyy = np.zeros(self.n * self.n)
        self.zzz = np.zeros(self.n * self.n)

        for i in range(0, self.n * self.n):
            self.xxx[i] = self.exes[i]
            self.yyy[i] = self.eyes[i]
            self.zzz[i] = self.ezes[i]
        self.places = np.vstack([self.xxx, self.yyy, self.zzz]).T

        self.Ui = RGI((self.X[:, 0, 0], self.Y[0, :, 0], self.Z[0, 0, :]), self.U, method='linear')
        self.Vi = RGI((self.X[:, 0, 0], self.Y[0, :, 0], self.Z[0, 0, :]), self.V, method='linear')
        self.Wi = RGI((self.X[:, 0, 0], self.Y[0, :, 0], self.Z[0, 0, :]), self.W, method='linear') 

    def simulate(self):
        self.ax = plt.figure().add_subplot(projection='3d')

        dt = self.step_size
        self.xlist = []
        self.ylist = []
        self.zlist = []

        for p in self.places:
            r = ode(self.IL3)
            r.set_integrator(self.integrmodel)
            lx = [p[0]]
            ly = [p[1]]
            lz = [p[2]]
            r.set_initial_value([lx[0], ly[0], lz[0]], 0)

            step = 0
            abort = False

            while r.successful():
                step += 1
                r.integrate(r.t + dt)
                x, y, z = r.y[0], r.y[1], r.y[2]

                lx.append(x)
                ly.append(y)
                lz.append(z)
                self.xlist.append(x)
                self.ylist.append(y)
                self.zlist.append(z)

                if self.a >= x or self.b <= x or self.a >= y or self.b <= y or self.a >= z or self.b <= z:
                    break

                if step >= self.max_step:
                    break

                if abort is True:
                    abort = False
                    break

            plt.plot(lx, ly, lz, 'k', lw=0.5)

        plt.show()

