import numpy as np
import math
import matplotlib.pyplot as plt


def solveimage(gridmin,gridlength,gridsize,beta,xgrad,ygrad,scale):

# def solveimage(beta,delxx,delyy,scale):
        redxpt=np.zeros(5)
        redypt=np.zeros(5)
        bluxpt=np.zeros(5)
        bluypt=np.zeros(5)
        image=np.zeros([21,2])


        x=np.linspace(gridmin[0],gridmin[0]+gridlength,gridsize)
        y=np.linspace(gridmin[1],gridmin[1]+gridlength,gridsize)

        gridstep=abs(x[2]-x[1])
#         print(gridstep)

        Y,X=np.meshgrid(x,y)
        delxx=-X+beta[0]+xgrad*scale
        delyy=-Y+beta[1]+ygrad*scale
        
#         plt.contour(X,Y,delyy,0,colors='blue')
#         plt.contour(X,Y,delxx,0,colors='red')
        
        im=0
        for i in range(1,len(x)):
            for j in range(1,len(y)):
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
                        
#                         print(xint,yint)

                        if xint>0.0 and xint<=1.0 and yint>0.0 and yint<=1.0:
                            im=im+1;
                            image[im-1,0], image[im-1,1]=gridmin[0]+((i-1)+xint)*gridstep, gridmin[1]+((j-1)+yint)*gridstep
#         print(image)
        img=image[~np.all(image == 0, axis=1)]
#         print(img)
        return img