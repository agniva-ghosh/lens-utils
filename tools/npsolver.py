import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.weave import converters

def solveimage(gridmin,gridlength,gridsize,beta,xgrad,ygrad,scale):

# def solveimage(beta,delxx,delyy,scale):
        redxpt=np.zeros(5)
        redypt=np.zeros(5)
        bluxpt=np.zeros(5)
        bluypt=np.zeros(5)
        image=np.zeros([21,2])


        x=np.linspace(gridmin,gridmin+gridlength,gridsize)
        y=np.linspace(gridmin,gridmin+gridlength,gridsize)

        gridstep=abs(x[2]-x[1])
#         print(gridstep)

        Y,X=np.meshgrid(x,y)
        delxx=-X+beta[0]+xgrad*scale
        delyy=-Y+beta[1]+ygrad*scale
        
#         plt.contour(X,Y,delyy,0,colors='blue')
#         plt.contour(X,Y,delxx,0,colors='red')





        
        code = '''
              double redA,redB,bluA,bluB,im[3,11],xint,yint;
              double bluxpt[2],bluypt[2],redxpt[2],redypt[2];
              int i,j,im,ix,iy,nminus,mminus;
              im=0;
              for(ix=2; ix<=1601; ix++)
              {
                for(iy=2; iy<=1601; iy++)
                {
                  nminus=0;

                  if( delxx[ix-1][iy-1]*delxx[ix][iy-1] < 0.0)
                  {
                    nminus=nminus+1;
                  }
                  if( delxx[ix][iy-1]*delxx[ix][iy] < 0.0)
                  {
                    nminus=nminus+1;
                  }
                  if( delxx[ix][iy]*delxx[ix-1][iy] < 0.0)
                  {
                    nminus=nminus+1;
                  }
                  if( delxx[ix-1][iy]*delxx[ix-1][iy-1] < 0.0)
                  {
                    nminus=nminus+1;
                  }
                  if(nminus==2)
                  {
                    mminus=0;
                    if( delyy[ix-1][iy-1]*delyy[ix][iy-1] < 0.0)
                    {
                      mminus=mminus+1;
                    }
                    if( delyy[ix][iy-1]*delyy[ix][iy] < 0.0)
                    {
                      mminus=mminus+1;
                    }
                    if( delyy[ix][iy]*delyy[ix-1][iy] < 0.0)
                    {
                      mminus=mminus+1;
                    }
                    if( delyy[ix-1][iy]*delyy[ix-1][iy-1] < 0.0)
                    {
                      mminus=mminus+1;
                    }
                    if(mminus==2)
                    {
                      i=0;
                      if(delxx[ix-1][iy-1]*delxx[ix][iy-1] <= 0.0)
                      {
                        i=i+1;
                        redxpt[i]=fabs(delxx[ix-1][iy-1]/(delxx[ix-1][iy-1]-delxx[ix][iy-1]));
                        redypt[i]=0.0;
                      }
                      if(delxx[ix][iy-1]*delxx[ix][iy] <= 0.0)
                      {
                        i=i+1;
                        redypt[i]=fabs(delxx[ix][iy-1]/(delxx[ix][iy-1]-delxx[ix][iy]));
                        redxpt[i]=1.0;
                      }
                      if(delxx[ix-1][iy]*delxx[ix][iy] <= 0.0)
                      {
                        i=i+1;
                        redxpt[i]=fabs(delxx[ix-1][iy]/(delxx[ix-1][iy]-delxx[ix][iy]));
                        redypt[i]=1.0;
                      }
                      if(delxx[ix-1][iy-1]*delxx[ix-1][iy] <= 0.0)
                      {
                        i=i+1;
                        redypt[i]=fabs(delxx[ix-1][iy-1]/(delxx[ix-1][iy-1]-delxx[ix-1][iy]));
                        redxpt[i]=0.0;
                      }

                      j=0;
                      if(delyy[ix-1][iy-1]*delyy[ix][iy-1] <= 0.0)
                      {
                        j=j+1;
                        bluxpt[j]=fabs(delyy[ix-1][iy-1]/(delyy[ix-1][iy-1]-delyy[ix][iy-1]));
                        bluypt[j]=0.0;
                      }
                      if(delyy[ix][iy-1]*delyy[ix][iy] <= 0.0)
                      {
                        j=j+1;
                        bluypt[j]=fabs(delyy[ix][iy-1]/(delyy[ix][iy-1]-delyy[ix][iy]));
                        bluxpt[j]=1.0;
                      }
                      if(delyy[ix-1][iy]*delyy[ix][iy] <= 0.0)
                      {
                        j=j+1;
                        bluxpt[j]=fabs(delyy[ix-1][iy]/(delyy[ix-1][iy]-delyy[ix][iy]));
                        bluypt[j]=1.0;
                      }
                      if(delyy[ix-1][iy-1]*delyy[ix-1][iy] <= 0.0)
                      {
                        j=j+1;
                        bluypt[j]=fabs(delyy[ix-1][iy-1]/(delyy[ix-1][iy-1]-delyy[ix-1][iy]));
                        bluxpt[j]=0.0;
                      }


                      redA=(redypt[1]-redypt[2])/(redxpt[1]-redxpt[2]);
                      redB= redypt[1]-redA*redxpt[1];
                      bluA=(bluypt[1]-bluypt[2])/(bluxpt[1]-bluxpt[2]);
                      bluB= bluypt[1]-bluA*bluxpt[1];

                      xint=(bluB-redB)/(redA-bluA);
                      yint=redA*xint+redB;


                      /*Verification of the True Image point*/

                      if(xint>0.0 && xint<=1.0 && yint>0.0 && yint<=1.0)
                      {
                        im=im+1;
                        im[1,im]=-gridmin+((ix-2)+xint)*gridstep;
                        im[2,im]=-gridmin+((iy-2)+yint)*gridstep;
                      }

                    }
                  }
                }
             }
             return_val = im;
        
        
            '''
        
        image = weave.inline(code,['delxx','delyy','gridstep','gridmin'],type_converters=converters.blitz,
                              compiler = 'gcc')
#         image[im-1,0], image[im-1,1]=gridmin[0]+((i-1)+xint)*gridstep, gridmin[1]+((j-1)+yint)*gridstep
#         print(image)
        image=np.array(image)
        img=image[~np.all(image == 0, axis=1)]
#         print(img)
        return img