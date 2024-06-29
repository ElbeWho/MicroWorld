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
            function: "free" or "wall"; optional argument informing how data wil be plotted,
                        with or without obstacle.

        ---------

        Velocity values prepared to be used at the end of the definition:
            u: values of velocity in the x direction
            v: values of velocity in the y direction
        """

        self.flag = function 
        self.arrows = []
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
            self.b = b
            self.mX, self.mY = np.mgrid[self.a:b:steps, -0.2:b:0.005]
            self.steps = self.mX.shape[0]
            self.u = 0.*self.mX
            self.v = 0.*self.mY 
            command = False
        if command == True:
            print("Possible boundary conditions are 'free' and 'wall'.")
        

    def stokeslet(self, r0: np.ndarray, F: np.ndarray, coef = 1):
        """
        Arguments:
            r0: position of the force
            F: direction and magnitude of the force
            coef: optional argument set to 1, can be used to multiply whole welocity field by choosen number
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
        self.arrows.append([r0[0], r0[1], F[0], F[1]])


    def dipole(self, r0: np.ndarray, F: np.ndarray, d: np.ndarray, coef = 1):
        """
        Arguments:
            r0: position of the force
            F: direction and magnitude of the force
            d: distance between forces
            coef: optional argument set to 1, can be used to multiply whole welocity field by choosen number
        """
        e = F/(F[0]**2+F[1]**2)**.5
        r=np.array([self.mX-r0[0], self.mY-r0[1]])
        modr=(r[0]**2+r[1]**2)**.5
        first_term = ((d[0]*r[0]+ d[1]*r[1])*e[:, np.newaxis, np.newaxis] - (e[0]*r[0]+e[1]*r[1])*d[:, np.newaxis, np.newaxis] )/modr**3    
        second_term =  -((d[0]*e[0]+d[1]*e[1])*r)/modr**3
        third_term = 3*((e[0]*r[0]+e[1]*r[1])*(d[0]*r[0]+ d[1]*r[1])*r)/modr**5 
        ud, vd = first_term  + second_term + third_term 
        self.u += ud*coef
        self.v += vd*coef
        if F[0]>=0 and F[1]>=0:
            self.arrows.append([r0[0], r0[1], F[0], F[1]])
            self.arrows.append([r0[0], r0[1], -F[0], -F[1]])
        else:
            self.arrows.append([-F[0], -F[1], F[0], F[1]])
            self.arrows.append([F[0], F[1], -F[0], -F[1]])

    def stokes_dipole(self, r0: np.ndarray, F: np.ndarray, coef = 1):
        """
        Arguments:
            r0: position of force
            F: direction and magnitude of the force
            coef: optional argument set to 1, can be used to multiply whole welocity field by choosen number
        """
        e = F/(F[0]**2+F[1]**2)**.5
        
        first_size = self.u.shape[0]
        second_size = self.u.shape[1]
        Id=np.zeros([second_size, first_size ])
        for i in range(0, second_size):
            Id[i,i]=1
        r=np.array([self.mX-r0[0], self.mY-r0[1]]) 
        modr=(r[0]**2+r[1]**2)**.5
        Idr = np.dot(r, Id) 
        second_term = (3*(e[0]*r[0]+e[1]*r[1])**2)*r/modr**5
        udip, vdip =    (-1)*Idr/modr**3  + second_term
        self.u += udip*coef
        self.v += vdip*coef
        self.arrows.append([r0[0], r0[1], F[0], F[1]])

    def rotlet(self, r0: np.ndarray, F: np.ndarray, d: np.ndarray, coef = 1):
        """
        Arguments:
            r0: position of the force
            F: direction and magnitude of the force
            d: distance between forces
            coef: optional argument set to 1, can be used to multiply whole welocity field by choosen number
        """
        e = F
        r = np.array([self.mX-r0[0], self.mY-r0[1]])
        modr = (r[0]**2+r[1]**2)**.5 
        ua, va = ((d[0]*r[0]+ d[1]*r[1])*e[:, np.newaxis, np.newaxis] - (e[0]*r[0]+e[1]*r[1])*d[:, np.newaxis, np.newaxis] )/modr**3
        self.u += ua*coef
        self.v += va*coef
        self.arrows.append([r0[0], r0[1], F[0], F[1]])

    def rotlet_R(self, r0: np.ndarray, R: np.ndarray, coef = 1):
        """
        Arguments:
            r0: position of the force
            R: torque acting at the origin
            coef: optional argument set to 1, can be used to multiply whole welocity field by choosen number
        """
        r = np.array([self.mX-r0[0], self.mY-r0[1]])
        modr = (r[0]**2+r[1]**2)**.5
        ua = -(R[2]*r[1])/modr**3 
        va = (R[2, np.newaxis, np.newaxis]*r[0])/modr**3
        self.u += ua*coef
        self.v += va*coef

    def stresslet(self, r0: np.ndarray, F: np.ndarray, d: np.ndarray, coef = 1):
        """
        Arguments:
            r0: position of the force
            F: direction and magnitude of the force
            d: distance between forces
            coef: optional argument set to 1, can be used to multiply whole welocity field by choosen number
        """
        e=F
        r=np.array([self.mX-r0[0], self.mY-r0[1]])
        modr=(r[0]**2+r[1]**2)**.5 
        dwa =  -((d[0]*e[0]+d[1]*e[1])*r)/modr**3
        trzy = 3*((e[0]*r[0]+e[1]*r[1])*(d[0]*r[0]+ d[1]*r[1])*r)/modr**5
        us, vs = dwa+trzy
        self.u += us*coef
        self.v += vs*coef
        self.arrows.append([r0[0], r0[1], F[0], F[1]])

    def source(self, r0: np.ndarray, P: float, coef=1):
        """
        Arguments:
            r0: position of the force
            M: magnitude of a source
            coef: optional argument set to 1, can be used to multiply whole welocity field by choosen number
        """
        r = np.array([self.mX-r0[0], self.mY-r0[1]])
        modr=(r[0]**2+r[1]**2)**.5 
        ur, vr  = P*r/modr**3
        self.u += ur*coef
        self.v += vr*coef


    def source_dipole(self, r0: np.ndarray, M: np.ndarray, coef = 1):
        """
        Arguments:
            r0: position of the force
            M: strenght of the dipole
            coef: optional argument set to 1, can be used to multiply whole welocity field by choosen number
        """
        r=np.array([self.mX-r0[0], self.mY-r0[1]])
        Id=np.array([[1,0],[0,1]])
        IdM=np.dot(Id, M) 
        rM=(r*M[:,np.newaxis,np.newaxis]).sum(axis=0)
        rrM=(r*rM[np.newaxis,])
        modr=(r[0]**2+r[1]**2)**.5
        second = 3*rrM/modr**5
        u0, v0 = -IdM[:,np.newaxis,np.newaxis]/modr**3+ second
        self.u += u0*coef
        self.v += v0*coef
        self.arrows.append([r0[0], r0[1], M[0], M[1]])

    def free_surf_par(self, r0: np.ndarray, Fpar: np.ndarray):
        """
        Arguments:
            r0: position of the force
            Fpar: direction and magnitude of the force, 
                    only parallel component will be considered
        """
        Frel = np.array([Fpar[1], 0])
        #real system
        self.stokeslet(r0, Frel)
        #image system
        self.stokeslet(-r0, Frel)
        self.arrows.append([r0[0], r0[1], Frel[0], Frel[1]])

    def free_surf_per(self, r0: np.ndarray, Fper: np.ndarray):
        """
        Arguments:
            r0: position of the force
            Fpar: direction and magnitude of the force, 
                    only perpendicular component will be considered
        """
        Frel = np.array([0, Fper[0]])
        #real system
        self.stokeslet(r0, Frel)
        #image system
        self.stokeslet(-r0, -Frel)
        self.arrows.append([r0[0], r0[1], Frel[0], Frel[1]])

    def hard_wall_par(self, r0: np.ndarray, Fpar: np.ndarray):
        """
        Arguments:
            r0: position of the force
            Fpar: direction and magnitude of the force, 
                    only perpendicular component will be considered
        """
        h = r0[1]
        rim =np.array([r0[0], -r0[1]])
        coef1 = 2*h
        coef2 = -2*h**2
        d = np.array([1,0])
        Frel = np.array([1, 0])*Fpar[0]
        Fadd = np.array([0, 1])
        #--real
        self.stokeslet(r0, Frel) 
        #--image
        self.stokeslet(rim, Frel, coef=-1 ) 
        self.dipole(rim, Fadd, d, coef=coef1)
        self.source_dipole(rim, Frel, coef=coef2)
        self.arrows.append([r0[0], r0[1], Frel[0], Frel[1]])

    def hard_wall_per(self, r0: np.ndarray, Fper: np.ndarray):
        """
        Arguments:
            r0: position of the force
            Fpar: direction and magnitude of the force, 
                    only perpendicular component will be considered
        """
        h = r0[1]
        rim =np.array([r0[0], -r0[1]])
        coef0 = -1
        coef1 = -2*h
        coef2 = 2*h**2
        d = np.array([0, 1])
        Frel = np.array([0, Fper[0]])
        #--real system
        self.stokeslet(r0, Frel) 
        #--image system
        self.stokeslet(rim, Frel, coef=coef0 ) 
        self.dipole(rim, Frel, d, coef=coef1)
        self.source_dipole(rim, Frel, coef=coef2)
    
    def totality_par(self, r0, Fpar: np.ndarray, ratio: float):
        """
        Arguments:
            r0: position of the force
            Fpar: direction and magnitude of the force, 
                    only perpendicular component will be considered
            ratio: ratio of fluids' viscosities
        """
        h = r0[1]
        rim =np.array([r0[0], -r0[1]])
        coef0 = (1-ratio)/(1+ratio)
        coef1 = (2*ratio*h)/(ratio+1)
        coef2 = -1*(2*ratio*h**2)/(ratio + 1)
        d = np.array([1,0])
        Frel = np.array([1, 0])*Fpar[0]
        Fadd = np.array([0, 1])
        #--real system
        self.stokeslet(r0, Frel) 
        #--image system
        self.stokeslet(rim, Frel, coef=coef0 ) 
        self.dipole(rim, Fadd, d, coef=coef1)
        self.source_dipole(rim, Frel, coef=coef2)
    
    def totality_per(self, r0, Fper, ratio):
        """
        Arguments:
            r0: position of the force
            Fpar: direction and magnitude of the force, 
                    only perpendicular component will be considered
            ratio: ratio of fluids' viscosities
        """
        h = r0[1]
        rim =np.array([r0[0], -r0[1]])
        coef0 = -1
        coef1 = -(2*ratio*h)/(ratio+1)
        coef2 = (2*ratio*h**2)/(ratio + 1)
        d = np.array([0, 1])
        Frel = np.array([0, Fper[0]])
        #--real system
        self.stokeslet(r0, Frel) 
        #--image system
        self.stokeslet(rim, Frel, coef=coef0 ) 
        self.dipole(rim, Frel, d, coef=coef1)
        self.source_dipole(rim, Frel, coef=coef2)

    def streamlines(self, xstart: np.ndarray, ystart: np.ndarray, mesh=False, direction=False, export=False, title='pass'):
        """
        Arguments:
            xstart: an array with starting points for x axis
            ystart: an array with starting points for y axis
        Optional:
            mesh: bool value, originally set to False, if set to True plots colormap 
            direction: bool value, originally set to False, if set to True shows directions 
                        and magnitude of forces with arrows
            export: bool value, originally set to False, if set to True function saves 
                        interpolated data of stremlines as .npz file
            title: sting value, originally set to "pass", if different string given plot 
                        will be saved with extension given within a string
        """
        plt.rcParams['text.usetex'] = True
        plt.rcParams.update({
            'font.size': 20,
            'text.usetex': True,
            'text.latex.preamble': r'\usepackage{dsfont}'
        })

        if self.flag == "free" and export==False:
            fig = plt.figure(figsize=(8, 8), facecolor="w")
            
        elif self.flag == "wall" and export==False:
            fig = plt.figure(figsize=(8, 4), facecolor="w")
        else:
            pass

        if mesh and export == False:
            ax = plt.axes()
            Z = np.sqrt(self.v**2 + self.u**2)
            self.image = ax.pcolormesh(self.mX.T, self.mY.T, Z.T,
                                        norm=colors.LogNorm(vmin=10**(-1), vmax=10**1),
                                        snap=True,
                                        cmap=plt.cm.inferno, rasterized=True,
                                        shading='gouraud', zorder=0)
            self.add_colorbar(self.image)
        elif mesh and export:
            print('Error: export does not support plotting. Mesh is set as True.')
        else:
            pass

        places = np.vstack([xstart, ystart]).T
        ui = RGI((self.mX[:, 0], self.mY[0, :]), self.u, method='linear')
        vi = RGI((self.mX[:, 0], self.mY[0, :]), self.v, method='linear')

        integrmodel = 'vode'
        max_step = 2000
        step_size = 0.01
        dt = step_size

        def Integrate_Line(t, coord):
            xi, yi = coord
            try:
                ex = ui([xi, yi])[0]
                ey = vi([xi, yi])[0]
            except:
                return [False, False]
            n = np.sqrt(ex**2 + ey**2)
            return [ex / n, ey / n]
        
        if export:
            xlist = []
            ylist = []
            zlist = []
            for p in places:
                r = ode(Integrate_Line)
                r.set_integrator(integrmodel)
                lx=[p[0]]
                ly=[p[1]]
                r.set_initial_value([lx[0], ly[0]], 0)
                step = 0
                while r.successful():
                    step += 1
                    r.integrate(r.t+dt)
                    x, y= r.y[0], r.y[1]
                    lx.append(x)
                    ly.append(y)
                    xlist.append(x)
                    ylist.append(y)
                    zlist.append(0)
                    if (x < self.a - 0.1 or x > self.b - 0.1 or y < self.a - 0.1 or y > self.b - 0.1 or
                    step >= max_step):
                        xlist.append(None)
                        ylist.append(None)
                        zlist.append(None)
                        break
            triple_list = []

            for i in range(0, len(xlist)):
                single_array = np.array([xlist[i], ylist[i], zlist[i]])
                print(single_array)
                triple_list.append(single_array)

            triple_array = np.array(triple_list)
            np.savez('output', triple_array)
        else:
            for p in places:
                r = ode(Integrate_Line)
                r.set_integrator(integrmodel)
                lx, ly = [p[0]], [p[1]]
                r.set_initial_value([lx[0], ly[0]], 0)
                step = 0
                while r.successful():
                    step += 1
                    r.integrate(r.t + dt)
                    x, y = r.y[0], r.y[1]
                    lx.append(x)
                    ly.append(y)
                    if (x < self.a - 0.1 or x > self.b - 0.1 or y < self.a - 0.1 or y > self.b - 0.1 or
                        step >= max_step):
                        break
                ax.plot(lx, ly, 'r' if step >= max_step else 'k')
                ax.set_xlim(self.a, self.b)
                ax.set_ylim(self.a, self.b)
                ax.set_aspect('equal')
                ax.set_xticks(np.arange(self.a, self.b + 1, 1))
                ax.set_yticks(np.arange(self.a, self.b + 1, 1))

        if direction and export==False:
            if len(self.arrows) == 0:
                print("For given singularity sollutions it is not possible to plot an arrow.")
            else:
                for i in range(0, len(self.arrows)):
                    origin = np.array([[self.arrows[i][0], self.arrows[i][1]]])
                    vector = np.array([[self.arrows[i][2], self.arrows[i][3]]])
                    print(self.arrows)
                    ax.quiver(origin[:, 0], origin[:, 1], vector[:, 0], vector[:, 1], scale=1,
                                angles='xy', scale_units='xy', color='b', zorder=5)
        elif direction and export:
            print('Eror: export mode does not support plotting. Direction is set as True.')
        else:
            pass
        if title != 'pass' and export==False:
            plt.savefig(title, bbox_inches='tight', pad_inches=0, dpi=200)
        elif title != 'pass' and export:
            print('Error: export does not support plotting. Title is given.')
        else:
            pass

    def add_colorbar(self, im, aspect=20, pad_fraction=0.5, **kwargs):
        """Adds a vertical color bar to an image plot."""
        current_ax = plt.gca()  # Get the current axes
        divider = axes_grid1.make_axes_locatable(current_ax)
        width = axes_grid1.axes_size.AxesY(current_ax, aspect=1. / aspect)
        pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
        cax = divider.append_axes("right", size=width, pad=pad)
        return im.axes.figure.colorbar(im, cax=cax, **kwargs)

    def plot(self, direction = False, title = 'pass'):
        """
        Optional:
            direction: bool value, originally set to False, if set to True shows directions 
                        and magnitude of forces with arrows
            title: sting value, originally set to "pass", if different string given plot 
                        will be saved with extension given within a string
        """
        
        plt.rcParams['text.usetex'] = True
        plt.rcParams.update({
                            'font.size': 20,
                            'text.usetex': True,
                            'text.latex.preamble': r'\usepackage{dsfont}'
                            })
        
        if self.flag == "free":
            fig = plt.figure(figsize=(8,8),facecolor="w")
            ax = plt.axes()
            ax.set_xlim(self.a, self.b)
            ax.set_ylim(self.a, self.b)
            ax.set_aspect('equal')
            ax.set_xticks(np.arange(self.a, self.b + 1, 1))
            ax.set_yticks(np.arange(self.a, self.b + 1, 1))
            
        elif self.flag == "wall":
            fig = plt.figure(figsize=(8,4),facecolor="w")
            ax = plt.axes()
        else:
            print("Error in boundary conditions.")

        Z = np.sqrt(self.v**2+self.u**2)
        
        self.image = ax.pcolormesh(self.mX.T, self.mY.T, Z.T,
                norm=colors.LogNorm(vmin= 10**(-3), vmax=10**1),
                #norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()),
                snap=True,
                cmap=plt.cm.inferno, rasterized=True, 
                shading='gouraud', zorder=0)
        
        plt.streamplot(self.mX.T, self.mY.T, self.u.T, self.v.T, 
               broken_streamlines=False, 
               density=0.3, 
               #z jakiegoś powodu nie działa dla rotlet 0.3
               color='k')

        if self.flag == "wall":
            plt.axhline(linewidth=10, y = -0.1, color=(0.5, 0.5, 0.5), linestyle = '-')
        else:
            pass
        
        self.add_colorbar(self.image)

        if direction== True:
            if len(self.arrows) == 0:
                print("For given singularity sollutions it is not possible to plot an arrow.")
            else:
                for i in range(0, len(self.arrows)):
                    origin = np.array([[self.arrows[i][0], self.arrows[i][1]]])
                    vector = np.array([[self.arrows[i][2], self.arrows[i][3]]])
                    print(self.arrows)
                    ax.quiver( origin[:, 0], origin[:, 1], vector[:, 0], vector[:, 1], scale=1,
                                angles='xy', scale_units='xy', color='b', zorder=5)
                
        else:
            pass

        if title != 'pass':
            plt.savefig(title, bbox_inches='tight', dpi=200)
        else:
            pass
    
    def show(self):
        """Shows image."""
        plt.show()
        