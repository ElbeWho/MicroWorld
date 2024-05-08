import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors 
import matplotlib.cbook as cbook
from matplotlib import cm
from mpl_toolkits import axes_grid1
from scipy.integrate import ode
from scipy.interpolate import RegularGridInterpolator as RGI

class Equations2D:
    """ 
    Stokes equation is relevant for the regime where Reynold's number Re << 1.
    This module provides ploting on meshgrid choosen singularity equations.
    It is possible to use build-in matplotlib function of streamlines or solve those by interpolating.
    """
    
    def __init__ (self, a: float, b: float, steps: float):
        """
        To initialize an object from this class, size of the plot must be specified.

        Arguments:
            a: lenght of the plot
            b: width of the plot
            steps: space between first and last point for both lenth and wisth

        ---------
            
        Velocity values prepared to be used at the end of the definition:
            u: values of velocity in the x direction
            v: values of velocity in the y direction
        """
        self.a = a
        self.b = b

        """xx = np.linspace(self.a, self.b, steps)
        yy = np.linspace(self.a, self.b, steps)
        self.mX, self.mY = np.meshgrid(xx,yy)
        self.u = np.zeros(self.mX.shape)
        self.v = np.zeros(self.mY.shape)"""
        self.mX, self.mY = np.mgrid[a:b:steps, a:b:steps]
        self.u = 0.*self.mX
        self.v = 0.*self.mY 
        #zeby dzialalo na mgrid trzeba jakos pozbyc sie zera
        

    def stokeslet(self, r0: np.ndarray, F: np.ndarray):
        """
        The fundamental sullution to the Stokes' equation with one point force applied
        is called Stokeslet.
        Point force defined as F = (8ηπr)δ(r)
        Equation:
        u(r)  = F(1 + r.r/|r|^2)

        Arguments:
            r0: position of the force
            F: direction and magnitude of the force
        """
        r=np.array([self.mX-r0[0], self.mY-r0[1]])
        Id=np.array([[1,0],[0,1]])
        Idf=np.dot(Id,F) 
        rTF=(r*F[:,np.newaxis,np.newaxis]).sum(axis=0)
        rrTF=(r*rTF[np.newaxis,])
        modr=(r[0]**2+r[1]**2)**.5
        u0, v0 =Idf[:,np.newaxis,np.newaxis]/modr+rrTF/modr**3.
        self.u += u0
        self.v += v0


    def dipole(self, r0: np.ndarray, F: np.ndarray, d: np.ndarray):
        """
        The first term of Tylor series for a Stokeslet at r with d beeing distance between two forces
        in opposite directions.
        Equation:
        u(r; F, d) = [(d x F) x r]/r^3 - [(d.F)r]/r^3 + [3(F.r)(d.r)r]/r^5
        
        Arguments:
            r0: position of the force
            F: direction and magnitude of the force
            d: distance between forces 
        """
        e = F/(F[0]**2+F[1]**2)**.5
        r=np.array([self.mX-r0[0], self.mY-r0[1]])
        modr=(r[0]**2+r[1]**2)**.5
        first_term = ((d[0]*r[0]+ d[1]*r[1])*F[:, np.newaxis, np.newaxis] - (e[0]*r[0]+e[1]*r[1])*d[:, np.newaxis, np.newaxis] )/modr**3    
        second_term =  -((d[0]*e[0]+d[1]*e[1])*r)/modr**3
        third_term = 3*((e[0]*r[0]+e[1]*r[1])*(d[0]*r[0]+ d[1]*r[1])*r)/modr**5 
        ud, vd = first_term  + second_term + third_term 
        self.u += ud
        self.v += vd

    def stokes_dipole(self, r0: np.ndarray, F: np.ndarray):
        """
        Dipole with stress tensor that is traceless and axisimetric. Both vector units of d and F are
        identical. This represents swimming flagellated bacteria.
        Equation:
        e = F/|F|
        u(r; S) = [-1/r^3 + (3(e.r)^2)/r^5]r
        Arguments:
            r0: position of force
            F: direction and magnitude of the force
        """
        e = F/(F[0]**2+F[1]**2)**.5
        Id=np.zeros([self.steps ,self.steps ])
        for i in range(0,self.steps ):
            Id[i,i]=1
        r=np.array([self.mX-r0[0], self.mY-r0[1]])
        modr=(r[0]**2+r[1]**2)**.5
        Idr = np.dot(r, Id)
        second_term = (3*(e[0]*r[0]+e[1]*r[1])**2)*r/modr**5
        udip, vdip =    (-1)*Idr/modr**3  + second_term
        self.u += udip
        self.v += vdip

    def rotlet(self, r0: np.ndarray, F: np.ndarray, d: np.ndarray):
        """
        Antysymmetric part of dipole velocity field.
        Equation:
        u(r; F, d) = [(d x F) x r]/r^3
        Arguments:
            r0: position of the force
            F: direction and magnitude of the force
            d: distance between forces 

        """
        e = F/(F[0]**2+F[1]**2)**.5
        r = np.array([self.mX-r0[0], self.mY-r0[1]])
        modr = (r[0]**2+r[1]**2)**.5 
        ua, va = ((d[0]*r[0]+ d[1]*r[1])*e[:, np.newaxis, np.newaxis] - (e[0]*r[0]+e[1]*r[1])*d[:, np.newaxis, np.newaxis] )/modr**3
        self.u += ua
        self.v += va
        
    def rotlet_R(self, r0: np.ndarray, R: np.ndarray):
        """
        Antysymmetric part of dipole velocity field with
        interpretation of the vector product R = d x F as torque. 
        Equation:
        u(r; R) = [R x r]/r^3
        Arguments:
            r0: position of the force
            R: torque acting at the origin
        """
        r = np.array([self.mX-r0[0], self.mY-r0[1]])
        modr = (r[0]**2+r[1]**2)**.5
        ua = -(R[2]*r[1])/modr**3 
        va = (R[2, np.newaxis, np.newaxis]*r[0])/modr**3
        self.u += ua
        self.v += va

    def stresslet(self, r0: np.ndarray, F: np.ndarray, d: np.ndarray,):
        """
        Symmetric part of dipole velocity field.
        Equation:
        u(r; F, d) = - [(d.F)r]/r^3 + [3(F.r)(d.r)r]/r^5
        Arguments:
            r0: position of the force
            F: direction and magnitude of the force
            d: distance between forces
        """
        e = F/(F[0]**2+F[1]**2)**.5
        r=np.array([self.mX-r0[0], self.mY-r0[1]])
        modr=(r[0]**2+r[1]**2)**.5 
        dwa =  -((d[0]*e[0]+d[1]*e[1])*r)/modr**3
        trzy = 3*((e[0]*r[0]+e[1]*r[1])*(d[0]*r[0]+ d[1]*r[1])*r)/modr**5 
        us, vs = dwa+trzy
        self.u += us
        self.v += vs

    def source(self, r0: np.ndarray, M: float):
        """
        Velocity field for a given logarytmic potental:
        φ(r) = M/(2π) ln(r),
        where M is the source magnitude.
        Equation:
        u(r) = Mr/(4πr^3)
        Arguments:
            r0: position of the force
            M: magnitude of a source
        """
        r = np.array([self.mX-r0[0], self.mY-r0[1]])
        modr=(r[0]**2+r[1]**2)**.5 
        ur, vr  = M*r/modr**3
        self.u += ur
        self.v += vr

    def source_dipole(self, r0: np.ndarray, M: np.ndarray):
        """
        Velocity field for a source dipole, also called source doublet.
        It is a source and a sink merging together.
        Equation:
        u(r; M) = 1/(4π) [1/r^3 + (3rr)/r^5] M
        Arguments:
            r0: position of the force
            M: strenght of the dipole
        """
        M = np.array([0,1])
        r=np.array([self.mX-r0[0], self.mY-r0[1]])
        Id=np.array([[1,0],[0,1]])
        IdM=np.dot(Id, M) 
        rM=(r*M[:,np.newaxis,np.newaxis]).sum(axis=0)
        rrM=(r*rM[np.newaxis,])
        modr=(r[0]**2+r[1]**2)**.5
        second = 3*rrM/modr**5
        u0, v0 = IdM[:,np.newaxis,np.newaxis]/modr**3- second
        self.u += u0
        self.v += v0

        ##--------------- tutaj będzie ścianka --------
        
    
    def streamlines(self):

        ui = RGI((self.mX[:, 0],self.mY[0, :]), self.u, method='linear')
        vi = RGI((self.mX[:, 0], self.mY[0, :]), self.v, method='linear')

        integrmodel = 'vode' # 'lsoda', 'dopri5', 'dop853'
        max_step = 2000
        step_size = 0.01
        fig, ax = plt.subplots()
        dt = step_size

        def Intergrate_Line(t, coord):
            # Function calculating Efield in every point for integrator
            global abort
            xi = coord[0]
            yi = coord[1]
            try:
                ex = ui([xi, yi])[0]
                ey = vi([xi, yi])[0]
                print(ui([xi, yi])[0])
            except:
                abort = True
                return [False, False]
            n = (ex**2+ey**2)**0.5
            return [ex/n, ey/n]

        xstart = np.linspace(-3.9, 3.9, 21)

        ystart = np.linspace(-3.9, -3.9, 21)
        places=np.vstack([xstart,ystart]).T

        for p in places:
            print("here")
            r = ode(Intergrate_Line)
            print("hree2")
            r.set_integrator(integrmodel)
            print("hree3")
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
                
                if self.a >= x or self.b <= x or self.a >= y or self.b <= y or \
                    step >= max_step:
                    break
                if step >= max_step:
                    
                    plt.plot(lx, ly, 'r') 
                else:
                    plt.plot(lx, ly, 'k')

        print("end now")
        

    def add_colorbar(self, im, aspect=20, pad_fraction=0.5, **kwargs):
        """Add a vertical color bar to an image plot."""
        current_ax = plt.gca()  # Get the current axes
        divider = axes_grid1.make_axes_locatable(current_ax)
        width = axes_grid1.axes_size.AxesY(current_ax, aspect=1. / aspect)
        pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
        cax = divider.append_axes("right", size=width, pad=pad)
        return im.axes.figure.colorbar(im, cax=cax, **kwargs)

    def __plot__(self):
        
        plt.rcParams['text.usetex'] = True
        plt.rcParams.update({
                            'font.size': 28,
                            'text.usetex': True,
                            'text.latex.preamble': r'\usepackage{dsfont}'
                            })
        fig = plt.figure(figsize=(8,8),facecolor="w")
        ax = plt.axes()
        Z = np.sqrt(self.v**2+self.u**2)
        
        
        self.image = ax.pcolormesh(self.mX, self.mY, Z,
                norm=colors.LogNorm(vmin= 10**(-2), vmax=10**1),
                #norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()),
                snap=True,
                cmap=plt.cm.inferno, rasterized=True, 
                shading='gouraud', zorder=0)


        plt.streamplot(self.mX.T, self.mY.T, self.u, self.v, 
               broken_streamlines=False, 
               density=0.35, 
               #z jakiegoś powodu nie działa dla rotlet 0.3
               color='k'
               )
        
        self.add_colorbar(self.image)

    
    def __show__(self):
        plt.show()
        
    def save_plot(self, title: str):
        """
        Arguments:
            title: title of the plot
        """
        plt.savefig(title + '.pdf', dpi=200)


