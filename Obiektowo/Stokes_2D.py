import numpy as np
import ast
import matplotlib.pyplot as plt
import matplotlib.colors as colors 
import matplotlib.cbook as cbook
from matplotlib import cm
import sympy
from mpl_toolkits import axes_grid1

class Equations2D:
    """ 
    Stokes equation is relevant for the regime where Reynold's number Re << 1.
    This module provides ploting on meshgrid choosen singularity equations. 
    """
    
    def __init__ (self, a: float, b: float, steps: float):
        """
        To initialize an object from this class size of the plot must be specified.

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
        self.steps = steps
        xx = np.linspace(-a, a, steps)
        yy = np.linspace(-b, b, steps)
        self.mX, self.mY = np.meshgrid(xx,yy)
        self.u = np.zeros(self.mX.shape)
        self.v = np.zeros(self.mY.shape)

    def stokeslet(self, r0: np.ndarray, f: np.ndarray):
        """
        The fundamental sullution to the Stokes' equation with one point force applied
        is called Stokeslet.
        Equation:
        V(r)  = F/(8ηπr) (1 + r.r/|r|^2)

        Arguments:
            r0: position of the force
            f: direction of the force
        """
    
        r=np.array([self.mX-r0[0], self.mY-r0[1]])
        Id=np.array([[1,0],[0,1]])
        Idf=np.dot(Id,f) 
        rTf=(r*f[:,np.newaxis,np.newaxis]).sum(axis=0)
        rrTf=(r*rTf[np.newaxis,])
        modr=(r[0]**2+r[1]**2)**.5
        u0, v0 =Idf[:,np.newaxis,np.newaxis]/modr+rrTf/modr**3.
        self.u += u0
        self.v += v0

    def dipole_Tylor(self, r0):

        Id=np.zeros([self.b ,self.b ])
        for i in range(0,self.b ):
            Id[i,i]=1
        
        r=np.array([self.mX-r0[0], self.mY-r0[1]])
    
        Idr = np.dot(r, Id)
        x2r = r*(r[1,:]**2)
        modr=(r[0]**2+r[1]**2)**.5
    
        udip, vdip =   3*x2r/modr**5  - Idr/modr**3
        self.u += udip
        self.v += vdip

    def dipole(self, r0, d, e):
        r=np.array([self.mX-r0[0], self.mY-r0[1]])
        modr=(r[0]**2+r[1]**2)**.5    
        dwa =  -((d[0]*e[0]+d[1]*e[1])*r)/modr**3
        trzy = 3*((e[0]*r[0]+e[1]*r[1])*(d[0]*r[0]+ d[1]*r[1])*r)/modr**5 
        jeden = ((d[0]*r[0]+ d[1]*r[1])*e[:, np.newaxis, np.newaxis] - (e[0]*r[0]+e[1]*r[1])*d[:, np.newaxis, np.newaxis] )/modr**3
        ud, vd = dwa+trzy+jeden
        self.u += ud
        self.v += vd
        

    def stresslet(self, r0, d, e):     
        r=np.array([self.mX-r0[0], self.mY-r0[1]])
        modr=(r[0]**2+r[1]**2)**.5    
        dwa =  -((d[0]*e[0]+d[1]*e[1])*r)/modr**3
        trzy = 3*((e[0]*r[0]+e[1]*r[1])*(d[0]*r[0]+ d[1]*r[1])*r)/modr**5 
        us, vs = dwa+trzy
        self.u += us
        self.v += vs

    def rotlet(self, r0, d, e):
        r = np.array([self.mX-r0[0], self.mY-r0[1]])
        modr = (r[0]**2+r[1]**2)**.5 
        jeden = ((d[0]*r[0]+ d[1]*r[1])*e[:, np.newaxis, np.newaxis] - (e[0]*r[0]+e[1]*r[1])*d[:, np.newaxis, np.newaxis] )/modr**3
        ua, va = jeden
        self.u += ua
        self.v += va

    def rotlet_R(self, r0, R):
        r = np.array([self.mX-r0[0], self.mY-r0[1]])
        modr = (r[0]**2+r[1]**2)**.5
        ua = -(R[2]*r[1])/modr**3 
        va = (R[2, np.newaxis, np.newaxis]*r[0])/modr**3
        self.u += ua
        self.v += va

    def source(self, r0):
        r = np.array([self.mX-r0[0], self.mY-r0[1]])
        modr=(r[0]**2+r[1]**2)**.5 
        macierz = r/modr**3
        ur, vr = macierz
        self.u += ur
        self.v += vr

    def source_doublet(self, e, r0):
        r = np.array([self.mX-r0[0], self.mY-r0[1]])
        modr=(r[0]**2+r[1]**2)**.5
        rer = 3*(r[0]*e[0]+r[1]*e[1])*r
        doublet = -rer/modr**5 + e[:,np.newaxis,np.newaxis]/modr**3
        usd, vsd = doublet
        self.u += usd
        self.v += vsd


    def source_dipole(self, r0):
        M = np.array([1,0])
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

    def Stokes_dipole(self, e, r0):
        Id=np.zeros([self.b ,self.b ])
        for i in range(0,self.b ):
            Id[i,i]=1
        
        r=np.array([self.mX-r0[0], self.mY-r0[1]])
        modr=(r[0]**2+r[1]**2)**.5

        Idr = np.dot(r, Id)

        second = (3*(e[0]*r[0]+e[1]*r[1])**2)*r/modr**5

        udip, vdip =   second - Idr/modr**3
        self.u += udip
        self.v += vdip
    

    def add_colorbar(self, im, aspect=20, pad_fraction=0.5, **kwargs):
        """Add a vertical color bar to an image plot."""
        current_ax = plt.gca()  # Get the current axes
        divider = axes_grid1.make_axes_locatable(current_ax)
        width = axes_grid1.axes_size.AxesY(current_ax, aspect=1. / aspect)
        pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
        cax = divider.append_axes("right", size=width, pad=pad)
        return im.axes.figure.colorbar(im, cax=cax, **kwargs)

    def __plot__(self):
        fig = plt.figure(figsize=(6,6),facecolor="w")
        ax = plt.axes()
        Z = np.sqrt(self.v**2+self.u**2)

        self.image = ax.pcolormesh(self.mX, self.mY, Z,
                #norm=colors.LogNorm(vmin= 10**(-1), vmax=10**1),
                norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()),
                snap=True,
                cmap=plt.cm.inferno, rasterized=True, 
                shading='gouraud', zorder=0)
        
        plt.streamplot(self.mX, self.mY, self.u, self.v, 
               broken_streamlines=False, 
               density=0.3, 
               #z jakiegoś powodu nie działa dla rotlet 0.3
               color='k'
               )
        
        self.add_colorbar(self.image)

       

        #self.image = ax.set_title(r'$\displaystyle\\v(r)='
         #      r'\frac{\mathbf{F}}{8 \pi \eta r } ( \mathds{1} + \frac{\mathbf{rr}}{r^2})$', fontsize=16, color='k')

    def __show__(self):
        plt.savefig('stokes_dipole_along_x.pdf', bbox_inches='tight', pad_inches=0, dpi=400)

        plt.show()
        #plt.savefig('fafik.pdf', bbox_inches='tight', pad_inches=0, dpi=400)
        
    #def save_plot(self):
     #    plt.savefig('whatttt.pdf', bbox_inches='tight', pad_inches=0, dpi=400)


