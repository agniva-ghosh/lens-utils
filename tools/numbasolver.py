from numba import jit,prange
import numpy as np
import math
from tools.lens import *

'''
Image Finding Algorithm using the lens equation residual input and linear interpolation
'''
@jit(nopython=True, parallel=True)
def solveimage(delxx,delyy,gridsize,image):
    redxpt=np.zeros(5, dtype=np.float64)
    redypt=np.zeros(5, dtype=np.float64)
    bluxpt=np.zeros(5, dtype=np.float64)
    bluypt=np.zeros(5, dtype=np.float64)

#         plt.contour(X,Y,delyy,0,colors='blue')
#         plt.contour(X,Y,delxx,0,colors='red')
#         @jit(nopython=True)
    im=0
    for i in range(1,gridsize):
        for j in range(1,gridsize):
            nminus=0

            if delxx[i-1,j-1]*delxx[i,j-1] < 0.0:
                nminus=nminus+1   

            if delxx[i,j-1]*delxx[i,j] < 0.0:   
                nminus=nminus+1   

            if delxx[i,j]*delxx[i-1,j] < 0.0:      
                nminus=nminus+1   

            if delxx[i-1,j]*delxx[i-1,j-1] < 0.0:
                nminus=nminus+1

            if nminus==2:
                mminus=0
                if delyy[i-1,j-1]*delyy[i,j-1] < 0.0:
                    mminus=mminus+1

                if delyy[i,j-1]*delyy[i,j] < 0.0:
                    mminus=mminus+1

                if delyy[i,j]*delyy[i-1,j] < 0.0:
                    mminus=mminus+1

                if delyy[i-1,j]*delyy[i-1,j-1] < 0.0:
                    mminus=mminus+1;

                if mminus==2:
                    k=0

                    if  delxx[i-1,j-1]*delxx[i,j-1] <= 0.0:
                        k+=1
                        redxpt[k]=abs(delxx[i-1,j-1]/(delxx[i-1,j-1]-delxx[i,j-1]))
                        redypt[k]=0.0

                    if  delxx[i,j-1]*delxx[i,j] <= 0.0:
                        k+=1
                        redypt[k]=abs(delxx[i,j-1]/(delxx[i,j-1]-delxx[i,j]))
                        redxpt[k]=1.0

                    if  delxx[i-1,j]*delxx[i,j] <= 0.0:
                        k+=1
                        redxpt[k]=abs(delxx[i-1,j]/(delxx[i-1,j]-delxx[i,j]))
                        redypt[k]=1.0

                    if  delxx[i-1,j-1]*delxx[i-1,j] <= 0.0:
                        k+=1
                        redypt[k]=abs(delxx[i-1,j-1]/(delxx[i-1,j-1]-delxx[i-1,j]))
                        redxpt[k]=0.0


                    l=0
                    if   delyy[i-1,j-1]*delyy[i,j-1] <= 0.0:
                        l+=1
                        bluxpt[l]=abs(delyy[i-1,j-1]/(delyy[i-1,j-1]-delyy[i,j-1]))
                        bluypt[l]=0.0

                    if  delyy[i,j-1]*delyy[i,j] <= 0.0:
                        l+=1
                        bluypt[l]=abs(delyy[i,j-1]/(delyy[i,j-1]-delyy[i,j]))
                        bluxpt[l]=1.0

                    if  delyy[i-1,j]*delyy[i,j] <= 0.0:
                        l+=1
                        bluxpt[l]=abs(delyy[i-1,j]/(delyy[i-1,j]-delyy[i,j]))
                        bluypt[l]=1.0

                    if  delyy[i-1,j-1]*delyy[i-1,j] <= 0.0:
                        l+=1
                        bluypt[l]=abs(delyy[i-1,j-1]/(delyy[i-1,j-1]-delyy[i-1,j]))
                        bluxpt[l]=0.0


                    redA=(redypt[1]-redypt[2])/(redxpt[1]-redxpt[2])
                    redB= redypt[1]-redA*redxpt[1]
                    bluA=(bluypt[1]-bluypt[2])/(bluxpt[1]-bluxpt[2])
                    bluB= bluypt[1]-bluA*bluxpt[1]

                    xint=(bluB-redB)/(redA-bluA);
                    yint=redA*xint+redB;


                    if xint>0.0 and xint<=1.0 and yint>0.0 and yint<=1.0:
                        
                        image[im,0]=(i-1)+xint
                        image[im,1]=(j-1)+yint
                        im=im+1
    return image


'''
function for finding the Einstein Radius. kappa=1 should mean it's the critical density??
'''
@jit(nopython=True)
def einsteinRadius(px,py,pr,kappa,area):
    
    pi=np.pi
    

    mass=np.zeros(len(pr), dtype=np.float64)
    a=np.zeros(len(pr), dtype=np.float64)

    for k in range(len(pr)):
        for i in range(len(px)):
            for j in range(len(py)):
                if px[i,j]**2+py[i,j]**2<=pr[k]**2: #check the flow! 
                    mass[k]+=kappa[i,j]*area
                    a[k]+=area

    for k in range(len(pr)):   
        if (mass[k]-pi*pr[k]*pr[k])**2==np.min((mass-pi*pr*pr)**2):
#             print(pr[k],mass[k]-pi*pr[k]*pr[k],np.min(abs(mass-pi*pr*pr)))
            thetaE=pr[k]

    return thetaE


'''
Calculation of Convergence, Shear and Deflection Angles maps from plummers
'''
@jit(nopython=True,parallel=True)
def mapPlummers(px,py,mass,pw,x,y,values):
    
    for i in range(len(x)):
        for j in prange(len(y)):
            
            values[0,i,j]=np.sum(4*mass*(pw*pw)/(((x[i]-px)*(x[i]-px)+(y[j]-py)*(y[j]-py)+pw*pw)**2)) #kappa

            values[1,i,j]=np.sum(4*mass*((y[j]-py)*(y[j]-py)-(x[i]-px)*(x[i]-px))/((x[i]-px)*(x[i]-px)+(y[j]-py)*(y[j]-py)+pw*pw)**2) #gamma1
            values[2,i,j]=np.sum(8*mass*((x[i]-px)*(y[j]-py))/((x[i]-px)*(x[i]-px)+(y[j]-py)*(y[j]-py)+pw*pw)**2) #gamma2

            values[3,i,j]=np.sum(4*mass*(x[i]-px)/((x[i]-px)*(x[i]-px)+(y[j]-py)*(y[j]-py)+pw*pw)) #gradx
            values[4,i,j]=np.sum(4*mass*(y[j]-py)/((x[i]-px)*(x[i]-px)+(y[j]-py)*(y[j]-py)+pw*pw)) #grady
            
    return values 


