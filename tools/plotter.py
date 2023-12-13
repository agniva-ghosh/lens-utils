import numpy as np
import math
import matplotlib.pyplot as plt


def createConvergenceMap(gridsize,kappa):
    
    gridlength=len(kappa)
    
    x=np.linspace(-gridsize,gridsize,gridlength)
    y=np.linspace(-gridsize,gridsize,gridlength)
    yy, xx = np.meshgrid(x,y)
    
    im = plt.pcolormesh(xx,yy,kappa)
    plt.contour(xx,yy,kappa,[1],colors='white')
    plt.colorbar(im)
    
    plt.gca().set_aspect('equal')
    plt.xlabel("x (arcsec)")
    plt.ylabel("y (arcsec)")
    
    
def createMagnificationMap(gridsize,mu):
    
    gridlength=len(mu)
    
    x=np.linspace(-gridsize,gridsize,gridlength)
    y=np.linspace(-gridsize,gridsize,gridlength)
    yy, xx = np.meshgrid(x,y)
    
    plt.contour(xx,yy,1/mu,1,colors='blue',linewidths=1,alpha=0.7)
    
    plt.gca().set_aspect('equal')
    plt.xlabel("x (arcsec)")
    plt.ylabel("y (arcsec)")

