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
    
    def __init__ (self, a: float, b: float, steps: float, function: str):
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
        
        self.flag = function 
        #słabe rozwiązanie z 'command' --> do przemyślenia

        command = True

        if self.flag == 'free':
            self.a = a
            self.b = b
            self.mX, self.mY = np.mgrid[a:b:steps, a:b:steps]
            self.steps = self.mX.shape[0]
            self.u = 0.*self.mX
            self.v = 0.*self.mY 
            print(self.steps)
            command = False
        
        if self.flag == 'wall':
            self.a = a - 0.4
            self.mX, self.mY = np.mgrid[self.a:b:steps, -0.2:b:0.005]
            print("mx: ", self.mX.shape)
            print("my: ", self.mY.shape)
            print("a: ", self.a)
            self.steps = self.mX.shape[0]
            self.u = 0.*self.mX
            self.v = 0.*self.mY 
            command = False
            
        if command == True:
            print("Possible functions are plot and export")
        

    def stokeslet(self, r0: np.ndarray, F: np.ndarray, coef = 1):
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
        self.u += u0*coef
        self.v += v0*coef

        #teraz to trzeba tak przerobic zeby stokeslet nie byl voidem


    def dipole(self, r0: np.ndarray, F: np.ndarray, d: np.ndarray, coef = 1):
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
        self.u += ud*coef
        self.v += vd*coef

    def stokes_dipole(self, r0: np.ndarray, F: np.ndarray, coef = 1):
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
        
        first_size = self.u.shape[0]
        second_size = self.u.shape[1]

        Id=np.zeros([second_size, first_size ])
        for i in range(0, second_size):
            Id[i,i]=1

        
        r=np.array([self.mX-r0[0], self.mY-r0[1]]) 



        #---the famous rest
        modr=(r[0]**2+r[1]**2)**.5
        Idr = np.dot(r, Id) 
        second_term = (3*(e[0]*r[0]+e[1]*r[1])**2)*r/modr**5
        udip, vdip =    (-1)*Idr/modr**3  + second_term
        self.u += udip*coef
        self.v += vdip*coef

    def rotlet(self, r0: np.ndarray, F: np.ndarray, d: np.ndarray, coef = 1):
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
        self.u += ua*coef
        self.v += va*coef
        
    def rotlet_R(self, r0: np.ndarray, R: np.ndarray, coef = 1):
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
        self.u += ua*coef
        self.v += va*coef

    def stresslet(self, r0: np.ndarray, F: np.ndarray, d: np.ndarray, coef = 1):
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
        self.u += us*coef
        self.v += vs*coef

    def source(self, r0: np.ndarray, M: float, coef=1):
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
        self.u += ur*coef
        self.v += vr*coef

    def source_dipole(self, r0: np.ndarray, M: np.ndarray, coef = 1):
        """
        Velocity field for a source dipole, also called source doublet.
        It is a source and a sink merging together.
        Equation:
        u(r; M) = 1/(4π) [1/r^3 + (3rr)/r^5] M
        Arguments:
            r0: position of the force
            M: strenght of the dipole
        """
        r=np.array([self.mX-r0[0], self.mY-r0[1]])
        Id=np.array([[1,0],[0,1]])
        IdM=np.dot(Id, M) 
        rM=(r*M[:,np.newaxis,np.newaxis]).sum(axis=0)
        rrM=(r*rM[np.newaxis,])
        modr=(r[0]**2+r[1]**2)**.5
        second = 3*rrM/modr**5
        u0, v0 = IdM[:,np.newaxis,np.newaxis]/modr**3- second
        self.u += u0*coef
        self.v += v0*coef

        ##--------------- tutaj będzie ścianka --------

    def free_suf_par(self, r0, Fpar, coef = 1):
        #to jest git
        Frel = np.array([Fpar[0], 0])
        erel = Frel/(Frel[0]**2+Frel[1]**2)*.5
        self.stokeslet(r0, erel, coef)
        eim = erel
        self.stokeslet(-r0, eim, coef)

    def free_suf_per(self, r0, Fper, coef=1):
        #to jest git tez
        Frel = np.array([0, Fper[0]])
        erel = Frel/(Frel[0]**2+Frel[1]**2)*.5
        self.stokeslet(r0, erel, coef)
        self.stokeslet(-r0, -erel, coef)

    def hard_wall_par(self, r0, Fpar ):

        ex = np.array([1, 0])
        ey = np.array([0, 1])

        h = r0[1]
        rim =np.array([r0[0], -r0[1]])
        coef1 = 2*h
        coef2 = -2*h**2
        d = np.array([1,0])
        radd = np.array([-r0[1], r0[0]])
        Frel = np.array([1, 0])
        Fadd = np.array([0, 1])
        
        #--real
        self.stokeslet(r0, Frel) 
        #--image
        self.stokeslet(rim, Frel, coef=-1 ) 
        self.dipole(rim, Fadd, d, coef=coef1)
        self.source_dipole(rim, -Frel, coef=coef2)

    def hard_wall_per(self, r0, Fper, ratio):
        versors = r0/(r0[0]**2+r0[1]**2)*.5
        ex = np.array([versors[0], 0])
        ey = np.array([0, versors[1]])
        h = r0[1] #as a scalar

        Frel = np.array([Fper[0], 0])
        #--real
        self.stokeslet(r0, Frel)
        #--image
        self.stokeslet(-r0, -Frel, coef=coef1)
        self.stokes_dipole(-r0, -Frel, coef=coef2)
        self.source_dipole(-r0, -Frel, coef=coef3)
    
    def totality(self, r0, Fper, ratio):

        ex = np.array([1, 0])
        ey = np.array([0, 1])
        h = r0[1]
        coef1 = (1-ratio)/(1+ratio)
        coef2 = (2*ratio*h)/(ratio+1)
        coef3 = -1*(2*ratio*h**2)/(ratio + 1)
        Frel = np.array([Fpar[0], 0])
        
        #--real
        self.stokeslet(ex, Frel)
        #--image
        self.stokeslet(-ex, -Frel, coef=coef1 )
        self.dipole(-ex, -Frel, -ey, coef=coef2)

        
    
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
        
        if self.flag == "free":

            fig = plt.figure(figsize=(6,6),facecolor="w")
        
        if self.flag == "wall":
            fig = plt.figure(figsize=(8,4),facecolor="w")

        ax = plt.axes()
        Z = np.sqrt(self.v**2+self.u**2)
        
        
        self.image = ax.pcolormesh(self.mX.T, self.mY.T, Z.T,
                norm=colors.LogNorm(vmin= 10**(-3), vmax=10**1),
                #norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()),
                snap=True,
                cmap=plt.cm.inferno, rasterized=True, 
                shading='gouraud', zorder=0)


        plt.streamplot(self.mX.T, self.mY.T, self.u.T, self.v.T, 
               broken_streamlines=False, 
               density=0.38, 
               #z jakiegoś powodu nie działa dla rotlet 0.3
               color='k'
               )
        
        if self.flag == "wall":
            plt.axhline(linewidth=8, y = -0.1, color=(0.5, 0.5, 0.5), linestyle = '-')
        
        self.add_colorbar(self.image)
        plt.savefig('wall.png', dpi=200)

    
    def __show__(self):
        plt.show()
        
    def save_plot(self, title: str):
        """
        Arguments:
            title: title of the plot
        """
        plt.savefig(title + '.pdf', dpi=200)


