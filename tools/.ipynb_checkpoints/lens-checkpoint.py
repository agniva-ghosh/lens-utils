import numpy as np
import math

class alphaPotLens:
    
    def gradx(self,x,y,b,alpha,s,q,K2): #Deflection Angle in x
        core = s*s+x*x+y*y/(q*q)+K2*x*y
        term=b*(alpha/2)*(2*x+K2*y)*core**(alpha/2-1)
        return term
        
    def grady(self,x,y,b,alpha,s,q,K2): #Deflection Angle in y
        core = s*s+x*x+y*y/(q*q)+K2*x*y
        term=b*(alpha/2)*(2*y/(q*q)+K2*x)*core**(alpha/2-1)   
        return term

    def gradxx(self,x,y,b,alpha,s,q,K2): 
        core = s*s+x*x+y*y/(q*q)+K2*x*y
        term1 = b*alpha*core**(alpha/2-1)
        term2 = b*alpha/2*(alpha/2-1)*((2*x+K2*y)**2)*core**(alpha/2-2)
        return term1+term2

    def gradyy(self,x,y,b,alpha,s,q,K2):
        core = s*s+x*x+y*y/(q*q)+K2*x*y
        term1 = b*alpha/(q*q)*core**(alpha/2-1)
        term2 = b*alpha/2*(alpha/2-1)*((2*y/(q*q)+K2*x)**2)*core**(alpha/2-2)
        return term1+term2
    
    def gradxy(self,x,y,b,alpha,s,q,K2):
        core = s*s+x*x+y*y/(q*q)+K2*x*y
        term1 = K2*b*alpha/2*core**(alpha/2-1)
        term2 = b*alpha/2*(alpha/2-1)*(2*x+K2*y)*(2*y/(q*q)+K2*x)*core**(alpha/2-2)
        return term1+term2
    
    def kappa(self,x,y,b,alpha,s,q,K2): #Convergence
        return (gradxx(x,y,b,alpha,s,q,K2)+gradyy(x,y,b,alpha,s,q,K2))/2
    
class plummerLens:
    
    def defx(self,x,y,pw,mass):
        core=(x*x+y*y)+pw*pw
        term=4*mass*x/core
        return term

    def defy(self,x,y,pw,mass):
        core=(x*x+y*y)+pw*pw
        term=4*mass*y/core
        return term
    
    def kappa(self,x,y,pw,mass):
        core=(x*x+y*y)+pw*pw
        term=4*mass*(pw*pw)*(core**(-2))
        return term
    
    def gamma1(self,x,y,pw,mass):
        core=(x*x+y*y)+pw*pw
        term=4*mass*(y*y-x*x)*(core**(-2)) #check the sign
        return term
        
    def gamma2(self,x,y,pw,mass):
        core=(x*x+y*y)+pw*pw
        term=8*mass*(x*y)*(core**(-2))
        return term
    
    
    
    
    