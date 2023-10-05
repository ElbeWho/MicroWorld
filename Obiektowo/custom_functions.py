import numpy as np
import ast
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.colors as colors 
import matplotlib.cbook as cbook
from matplotlib import cm
import sympy
from mpl_toolkits import axes_grid1

class Stokes:
    
    def __init__(self, a, b ):
        a=5
        b=100
        self.a=a
        self.b= b
        
        xx = np.linspace(-a, a, b)
        yy = np.linspace(-a, a, b)

        mX, mY = np.meshgrid(xx,yy)
        
        plt.rcParams['text.usetex'] = True
        plt.rcParams.update({
        'font.size': 8,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{dsfont}'
            })
        
#czy ja moge chciec miec same dane bez plotowania ich na wykresie
#czy robie z tego stricte pod plotting
#ale wsm to taki jest cel robic plots
#mX, mY lecą w konstruktorze więc jest bezpiecznie - wcale nie jest bo ich nie widac xd but whyyyyyy


    def stokeslet(f, r0):
        Id=np.array([[1,0],[0,1]])
        r=np.array([mX-r0[0],mY-r0[1]])

        Idf=np.dot(Id,f) 
    
        rTf=(r*f[:,np.newaxis,np.newaxis]).sum(axis=0)
        rrTf=(r*rTf[np.newaxis,])
        modr=(r[0]**2+r[1]**2)**.5
    
        u,v=Idf[:,np.newaxis,np.newaxis]/modr[np.newaxis]+rrTf/modr**3.
        
        return [u,v]
    
    def B_dir(t,p,fx,fz):
        ex = fx(p[0],p[1])
        ez = fz(p[0],p[1])
        n = (ex**2.0+ez**2.0)**0.5
        return [ex/n, ez/n]
    
    def add_colorbar(im, aspect=20, pad_fraction=0.5, **kwargs):
        
        divider = axes_grid1.make_axes_locatable(im.axes)
        width = axes_grid1.axes_size.AxesY(im.axes, aspect=1./aspect)
        pad = axes_grid1.axes_size.Fraction(pad_fraction, width)
        current_ax = plt.gca()
        cax = divider.append_axes("right", size=width, pad=pad)
        plt.sca(current_ax)
        return im.axes.figure.colorbar(im, cax=cax,**kwargs)

r0=np.array([0,0])            # position of the force red arrow - touchdown point
f=np.array([0,1]) 


