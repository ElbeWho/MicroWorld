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
    def __init__ (self, plot_title):
        self.plot_title = plot_title

    def get_name(self):
        print(self.plot_title)

    def font(self): #nie wiem czy będę umiała tego uyc w przyszłości
        plt.rcParams['text.usetex'] = True
        plt.rcParams.update({
        'font.size': 8,
        'text.usetex': True,
        'text.latex.preamble': r'\usepackage{dsfont}'
        })

    def __meshgrid__(self, a, b):
        self.a = a
        self.b = b
        xx = np.linspace(-a, a, b)
        yy = np.linspace(-a, a, b)
        self.mX, self.mY = np.meshgrid(xx,yy)
    
    def entries(self, entries):
        self.entries = entries
    
    def get_entries(self):
        print(self.entries)
  

    def stokeslet(self, f,r0):
        Id=np.array([[1,0],[0,1]])
        r=np.array([self.mX-r0[0], self.mY-r0[1]])

        Idf=np.dot(Id,f) #mnożenie macierzy f oraz Id
        
        rTf=(r*f[:,np.newaxis,np.newaxis]).sum(axis=0)
    
        rrTf=(r*rTf[np.newaxis,])
        modr=(r[0]**2+r[1]**2)**.5
    
        u, v =Idf[:,np.newaxis,np.newaxis]/modr+rrTf/modr**3.
        self.u = u
        self.v = v

    def many_stokeslets(self):
        for j in range(0, len(self.entries)):
            self.stokeslet( self.entries[j][0], self.entries[j][1])
            self.u += self.u
            self.v += self.v

    
    def get_velocities(self):
        print("v: ", self.v)
        print("u: ", self.u)

    def __plot__(self):
        fig = plt.figure(figsize=(6,6),facecolor="w")
        ax = plt.axes()
        Z = np.sqrt(self.v**2+self.u**2)

        self.image = ax.pcolormesh(self.mX, self.mY, Z,
                norm=colors.LogNorm(vmin= 10**(-1), vmax=10**1),
                #norm=colors.LogNorm(vmin=Z.min(), vmax=Z.max()),
                snap=True,
                cmap=plt.cm.inferno, rasterized=True, 
                shading='gouraud', zorder=0)
        
        
    def save_plot(self):
        plt.savefig('monopole_title.pdf', bbox_inches='tight', pad_inches=0, dpi=400)
        plt.show()