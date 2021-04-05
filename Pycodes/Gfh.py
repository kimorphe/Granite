import numpy as np
import matplotlib.pyplot as plt


class Gint:
    def __init__(self,cp):
        self.cp=cp
    def eval(self,f,h,x):
        r=np.sqrt(x*x+h*h)
        ka=2.*np.pi*f/self.cp;
        fx=np.exp(1j*r*ka)/np.sqrt(r)
        return(fx)
    def eval_reg(self,f,x):
        ka=2.*np.pi*f/self.cp;
        fx=(np.exp(1j*x*ka)-1.)/np.sqrt(x+1.e-20)
        return(fx)
    def intg(self,f,h,w,npnt):
        x=np.linspace(0,w*0.5,npnt)
        dx=x[1]-x[0];
        fx=self.eval(f,h,x)
        fx[0]*=0.5
        fx[-1]*=0.5
        J=np.sum(fx)*dx
        return(J)

    def sintg(self,f,w,npnt):
        x=np.linspace(0,w*0.5,npnt)
        dx=x[1]-x[0];
        fx=self.eval_reg(f,x)
        fx[0]*=0.5
        fx[-1]*=0.5
        J=np.sum(fx)*dx
        J+=np.sqrt(2*w)
        return(J)
    def Decay(self,freq,cps,h,w,npnt):
        Dcy=[]
        k=0
        for f in freq:
            self.cp=cps[k]
            J=self.intg(f,h,w,npnt)
            J0=self.sintg(f,w,npnt)
            Dcy.append(np.abs(J/J0))
            k+=1
        return(np.array(Dcy))

if __name__=="__main__":

    w=2.0   # [mm]
    cp=5.0  # [km/s]
    f=0.5   # [MHz] 
    Np=301;

    gw=Gint(cp)
    hs=np.linspace(0.001,3.5,41)   # [mm]
    J=[]
    for h in hs:
        J.append(gw.intg(f,h,w,Np))
    J0=gw.sintg(f,w,Np)

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(hs,np.real(J),"-r",label="real")
    ax.plot(hs,np.imag(J),"-b",label="imag")
    ax.plot(hs,np.abs(J),"-k",linewidth=2,label="magnitude")
    ax.legend()

    ax.plot(0,np.real(J0),"ro")
    ax.plot(0,np.imag(J0),"bo")
    ax.plot(0,np.abs(J0),"ko")
    ax.grid(True)
    ax.set_xlabel("h [mm]")
    ax.set_ylabel("u(0,h)")

    freq=np.linspace(0.1,3.0,31)

    Dcy=[]
    Dcy2=[]
    h0=3.42
    w=2.0   # [mm]
    cp=5.0  # [km/s]
    for f in freq:
        gw.cp=cp
        uz=gw.intg(f,h0,w,Np)
        u0=gw.sintg(f,w,Np)
        Dcy.append(abs(uz/u0))
    
    cps=np.ones(len(freq))*cp
    Dcy2=gw.Decay(freq,cps,h0,w,Np)

    fig2=plt.figure()
    bx=fig2.add_subplot(111)
    bx.plot(freq,Dcy,"-ko",markersize=10)
    bx.plot(freq,Dcy2,"-g")
    bx.grid(True)
    bx.set_xlabel("frequency [MHz]")
    bx.set_ylabel("decay ratio |u(0,h)/u(0,0)|")
    bx.set_ylim([0,0.6])


    plt.show()

