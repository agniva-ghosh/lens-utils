class Irtysh:

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