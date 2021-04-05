# -*- coding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import data_dir as ddr 
import Gfh

def sigmoid(x,x90):
    a=np.log(0.9/(1-0.9))/x90
    return(1/(1+np.exp(-a*x)))

class AWV:
    def __init__(self):
        self.cp_ready=False
        self.s2_ready=False
    def blank(self,ndat):
        self.Nt=ndat
        self.t0=0.0
        self.amp=np.zeros(ndat)
        self.time=np.arange(ndat)
        self.dt=1
    def burst(self,nbrst,f0,Nt,Td):
        self.Nt=Nt
        self.time=np.linspace(0,Td,Nt)
        self.dt=self.time[1]-self.time[0]
        self.amp=np.zeros(self.Nt)
        i0=int(1/f0/self.dt*nbrst)
        self.amp[0:i0]=np.cos(2*np.pi*f0*self.time[0:i0])
        df=1/(self.dt*Nt)
        self.freq=df*np.arange(Nt)
        self.freq=df

    def Amod_Gauss(self,t0,sig):
        arg=(self.time-t0)/sig;
        arg*=arg
        arg*=0.5
        Env=np.exp(-arg)
        self.amp*=Env

    def set_taxis(self,t1,t2,ndat=0):
        if ndat != 0:
            Nt=ndat
            self.Nt=Nt
        self.time=np.linspace(t1,t2,Nt)
        self.dt=time[1]-time[0]
        self.df=1/(dt*Nt)
        self.freq=df*np.arange(Nt)
    def get_fnum(self,f0):
        return(np.argmin(abs(self.freq-f0)))
    def load(self,fname,stashed=False):
        self.cp_ready=False
        self.s2_ready=False
        fp=open(fname,"r")
        time=[]
        amp=[]
        for row in fp:
            dat=row.strip().split(",")
            time.append(float(dat[0]))
            amp.append(float(dat[1]))
        amp=amp-np.mean(amp)
        amp=amp*2.0
        time=np.array(time)*1.e6
        amp=np.array(amp)
        #amp-=np.mean(amp)
        dt=time[1]-time[0];
        Nt=len(time);
        df=1/(dt*Nt)
        freq=df*np.arange(Nt)
        fp.close()

        self.amp=amp
        self.time=time
        self.freq=freq
        self.dt=dt
        self.df=df
        self.Nt=len(amp)

        if stashed:
            self.amp_raw=np.ones(self.Nt)
            self.amp_raw+=self.amp

    def FFT(self):
        Amp=np.fft.fft(self.amp);
        self.Amp=Amp
    def Hilbert(self,FFT=True):
        if FFT:
            self.Amp=np.fft.fft(self.amp)
        nt=self.Nt
        nt2=int(nt/2)
        Amp=np.zeros(nt).astype(complex)
        Amp+=self.Amp
        Amp[nt2:nt]*=-1j
        amph=np.real(np.fft.ifft(Amp))
        return(amph)

    def Butterworth(self,tb,w_6dB,mexp,apply=True):
        tt=self.time-tb;
        Phi=1+(tt/w_6dB)**mexp
        self.Phi=1/Phi
        if apply:
            self.amp*=self.Phi
        return(self.Phi)
    def Gauss(self,tb,sig,apply=True):
        arg=(self.time-tb)/sig
        arg*=arg
        self.Phi=np.exp(-arg*0.5)
        if apply:
            self.amp*=self.Phi
        return(self.Phi)
    def plot_ascan(self,ax,name=""):
        ax.plot(self.time,self.amp,label=name)
        ax.grid(True)
    def plot_FFT(self,ax,FFT=True,name=""):
        if FFT:
            self.Amp=np.fft.fft(self.amp)
        ax.semilogy(self.freq,abs(self.Amp),label=name)
        ax.set_xlim([self.freq[1],3])
        ax.set_ylim([1.e-03,10])
        ax.grid(True)
    def Wiener(self,eps=0.2,Apply=True,FFT=True):
        if FFT:
            self.Amp=np.fft.fft(self.amp)
        Smax=np.max(np.abs(self.Amp))
        S2=abs(self.Amp)/Smax
        S2*=S2
        Wnr=S2/(S2+eps*eps)
        if Apply:
            self.Amp*=Wnr
        return(Wnr)
    def gdelay(self,FFT=True):
        if FFT:
            self.Amp=np.fft.fft(self.amp)
        Nt=self.Nt
        phi=-np.unwrap(np.angle(self.Amp))
        tg=np.zeros(Nt)
        tmp=np.diff(phi)
        tg[0:Nt-1]+=tmp;
        tg[1:Nt]+=tmp;
        tg[1:Nt-1]/=2.0;
        tg/=(2.*np.pi*self.df)
        self.tg=tg

    def plot_gdelay(self,ax,FFT=True,name=""):
        if FFT:
            self.Amp=np.fft.fft(self.amp)
        Nt=self.Nt
        phi=-np.unwrap(np.angle(self.Amp))
        tg=np.zeros(Nt)
        tmp=np.diff(phi)
        tg[0:Nt-1]+=tmp;
        tg[1:Nt]+=tmp;
        tg[1:Nt-1]/=2.0;
        tg/=(2.*np.pi*self.df)
        self.tg=tg

        ax.plot(self.freq,tg,".",markersize=3,label=name)
        ax.set_xlim([0,3])
        ax.grid(True)
    def decay_factor(self,H,FFT=True,bT=1.0):
        if FFT:
            self.Amp=np.fft.fft(self.amp)
        k2=-np.log(np.abs(self.Amp*bT))/H
        self.s2=k2/(self.freq*2.*np.pi)
        self.s2_ready=True
    def phase_vel(self,H,FFT=True):
        if FFT:
            self.Amp=np.fft.fft(self.amp)
        f1=0.5; 
        nf1=self.get_fnum(f1)

        phi_h=-np.unwrap(np.angle(self.Amp[nf1:self.Nt+1]))
        phi_l=-np.unwrap(np.angle(self.Amp[0:nf1+1]))
        phi_l=np.linspace(0.00001,phi_h[0],nf1+1)
        Phi=np.hstack([phi_l[0:nf1],phi_h]);
        omg=self.freq*np.pi*2
        cp=H*omg/Phi;
        self.cp=cp
        self.cp_ready=True
    def Et(self,rho=1.0):
        if not self.cp_ready:
            print("phase velocity not available")
            return()
        if not self.s2_ready:
            print("s2 (decay factor) not available")
            return()
        s=1./self.cp+1j*self.s2
        s=1./self.cp
        s=1j*self.s2
        z2=1/(s*s)
        z2[0]=0
        Nt=self.Nt
        Nt2=int(Nt/2)
        z2[Nt2:Nt+1]=0
        z2*=2;
        z2*=self.Wiener(eps=0.1,FFT=False,Apply=False)
        dt=self.time[1]-self.time[0]
        df=self.freq[1]-self.freq[0]
        self.Gw=z2*rho
        Gt=2*np.fft.ifft(z2)*rho*self.Nt*df
        self.Gt=Gt
        #self.Gt[Nt2:Nt+1]=0
        tmp=np.cumsum(np.real(Gt))*dt
        Ft=np.zeros(self.Nt)
        Ft[1:self.Nt+1]=tmp[0:-1]
        return(Ft)

    def WinParam(self):
        t0=11; x90=1
        Zmd=sigmoid(self.time-t0,x90)
        self.amp*=Zmd
        tb=12   # assumed mean TOF
        Sig2=6  # 2-sigma width
        Sig=Sig2/2
        self.Gauss(tb,Sig,apply=True)   # Preliminary windowing
        amph=self.Hilbert(FFT=True) # Hilber transform
        zamp=abs(self.amp+1j*amph)  # analytic signal
        St=np.cumsum(abs(zamp))     # cummulative envelope
        St/=np.max(St)              # normalization
        t50=self.time[ np.argmin(abs(St-0.5))]  # get t10
        t10=self.time[ np.argmin(abs(St-0.1))]  # get t50
        t90=self.time[ np.argmin(abs(St-0.9))]  # get t90
        tb=t50; # Window parameter (mean)
        sig=t50-t10 # Window parameter (std)
        sig*=1.0
        return(tb,sig,t90)

if __name__=="__main__":


    #   ---- Reference Waveform ------
    aref=AWV()  # reference wave container
    dref=ddr.DATA_DIR("Ref")
    dref.show()
    tb_ref,tw_ref,mexp=dref.get_win_params()
    fnref="1MHznew.csv"
    aref.load(fnref)

    ht=3.42 # [mm] sample thickness
    gw=Gfh.Gint(5.0)
    wd=2    # [mm] source width (wedge tip width)
    awv=AWV()   # A-scan wave container

    Mineral="Na" # Mineral type 
    Mineral="Qt" # Mineral type 
    Mineral="K" # Mineral type 

    #   ---- Transmitted Waveform -------
    ddir=ddr.DATA_DIR(Mineral)
    ddir.show()
    dir_name=ddir.get_dir_name()
    tb,tw_6dB,mexp=ddir.get_win_params()
    nfile=ddir.get_nfile()
    bT=ddir.get_bT()

    num=201     # File No.
    num=345
    num=351
    num=203     # File No.
    fname=dir_name+"/scope_"+str(num)+".csv"
    awv.load(fname) # load A-scan

    #   --- A-scan Presentations ---
    fig1=plt.figure()
    ax=fig1.add_subplot(111)
    aref.plot_ascan(ax,name=fnref)
    awv.plot_ascan(ax,name=fname)
    #awv.Butterworth(tb,tw_6dB,mexp,apply=True)  # Apply Butterworth window
    #aref.Butterworth(tb_ref,tw_ref,mexp,apply=True) # Apply Butterworth window
    (t50,sig,t90)=aref.WinParam()
    aref.Gauss(t50,sig,apply=True)
    (t50,sig,t90)=awv.WinParam()
    awv.Gauss(t50,sig,apply=True)
    aref.plot_ascan(ax,name=fnref+"(windowed)")
    awv.plot_ascan(ax,name=fname+"(windowed)")

    #   --- Fourier Amplitude & Group Delay ----
    fig2=plt.figure()
    bx=fig2.add_subplot(211)
    cx=fig2.add_subplot(212)

    aref.plot_FFT(bx,FFT=True,name=fnref)  # FFT + plotting
    awv.plot_FFT(bx,FFT=True,name=fname)   # FFT + plotting
    Wnr=awv.Wiener(eps=0.2,FFT=False)
    awv.Amp/=aref.Amp   # Deconvolution
    awv.plot_FFT(bx,FFT=False,name="transfer func.")   # FFT + plotting
    bx.semilogy(awv.freq,Wnr,label="Winer")

    awv.plot_gdelay(cx,FFT=False,name="transfer func.")
    aref.plot_gdelay(cx,FFT=False,name=fnref)
    
    awv.phase_vel(ht,FFT=False) #,Gdelay=False)
    Gm=gw.Decay(awv.freq,awv.cp,ht,wd,301)
    nf1=awv.get_fnum(0.5)
    nf2=awv.get_fnum(2.0)
    print("Gm^-1=",np.mean(1/Gm[nf1:nf2]))
    #awv.Amp=awv.Amp/Gm

    Z1=dref.rho*dref.cp
    Z2=ddir.rho*awv.cp
    Tcoef=2*Z1/(Z1+Z2)
    B=1.77
    #awv.Amp=awv.Amp/Tcoef/B
    awv.Amp[0]=0.0
    awv.plot_FFT(bx,FFT=False,name="Gw compensated")   # plotting

    fig3=plt.figure()
    dx1=fig3.add_subplot(211)
    dx2=fig3.add_subplot(212)

    dx1.plot(awv.freq,ht/awv.tg,label="gourp velocity")
    dx1.set_ylim([0,8])
    dx1.plot(awv.freq,awv.cp,"-",label="phase velocity")
    awv.decay_factor(ht,FFT=False,bT=1.0)
    #dx2.plot(awv.freq,awv.s2*awv.cp,".-")
    #dx2.plot(awv.freq,awv.s2,".-")

    fig4=plt.figure()
    ex=fig4.add_subplot(111)
    phi=2*np.arctan(awv.s2*awv.cp)
    s1=1/awv.cp; s1[0]=1
    s2=awv.s2; s2[0]=1
    indx=np.argwhere(s2<0)
    print(indx)
    s=s1+1j*s2
    z2=1/(s*s)
    z2*=awv.Wiener(eps=0.1,FFT=True,Apply=False)
    """
    dt=awv.time[1]-awv.time[0]
    df=awv.freq[1]-awv.freq[0]
    Gt=2*np.fft.ifft(z2)*ddir.rho*awv.Nt*df
    Ft=np.cumsum(np.real(Gt))*dt
    """
    Ft=awv.Et(rho=ddir.rho)
    ex.plot(awv.time,np.real(awv.Gt))
    ex.plot(awv.time,Ft)
    #dx2.plot(awv.freq,np.abs(z2))

    uin=AWV()
    f0=1.0; Nt=awv.Nt;Td=awv.time[-1];
    nbrst=4
    uin.burst(nbrst,f0,Nt,Td)
    T0=1/f0; tb=T0*nbrst*0.5; nsig=6;sig=nbrst*T0/nsig
    uin.Amod_Gauss(tb,sig)
    #uin.FFT()
    uin.Amp=np.fft.ifft(uin.amp)
    nt2=int(uin.Nt*0.5)
    uin.Amp[nt2:uin.Nt+1]=0
    indx=np.argwhere(awv.s2<0)
    s2=awv.s2
    s2[indx]=0
    #s2[:]=0
    omg=2*np.pi*awv.freq
    nf1=awv.get_fnum(0.5)
    kin=(1/awv.cp+1j*s2)*omg-1j*omg[nf1]*s2[nf1]
    kin[0]=0
    xcod=[10,20,40,60,80,100]
    ofst=0

    f0=1.0
    fsig=0.4
    arg=(awv.freq-f0)/fsig
    W=np.exp(-arg*arg*0.5)
    for xx in xcod:
        #Uw=np.exp(1j*kin*xx)*np.conj(uin.Amp)*2
        #Uw=np.exp(1j*kin*xx)*(uin.Amp)*2
        Uw=np.exp(1j*kin*xx)*(W)*2
        #ut=np.fft.fft(Uw)/uin.Nt
        ut=np.fft.fft(Uw)#/uin.Nt
        dx2.plot(uin.time,np.real(ut)+ofst)
        #ofst+=20
    #dx2.set_xlim([0.2,3])


    fig5=plt.figure()
    ax_fig5=fig5.add_subplot(111)
    uin=np.real(np.fft.ifft(W))
    nt2=int(aref.Nt/2)
    T2=awv.time[nt2]
    nt2d=awv.Nt-nt2
    indx=range(nt2)
    jndx=np.arange(nt2,awv.Nt,1)
    I=np.hstack([jndx,indx])
    U=np.zeros(awv.Nt)

    U[0:nt2d]+=uin[jndx]
    U[nt2d:awv.Nt+1]+=uin[indx]
    #ax_fig5.plot(aref.time-T2,uin)
    #ax_fig5.plot(aref.time-T2,uin[I]/np.max(uin),label="in")
    ax_fig5.grid(True)

    Ew=awv.Gw*W*2
    Et=np.fft.ifft(Ew)
    #ax_fig5.plot(aref.time-T2,np.real(Et[I])/np.max(Et),label="out")
    #ax_fig5.legend()
    ax_fig5.plot(awv.freq,np.abs(awv.Gw))
    #ax_fig5.plot(awv.freq,np.degrees(np.angle(awv.Gw)))
    ax_fig5.set_xlim([0,3])

    #-------------- Graph appearances -----------------
    fsz=14
    axes=[ax,bx,cx,dx1,dx2,ex]
    for axis in axes:
        axis.tick_params(labelsize=fsz)
        axis.grid(True)

    dx1.set_xlim([0.4,2.0])
    #dx2.set_ylim([3,10])
    dx2.set_xlabel("frequency [MHz]",fontsize=fsz)
    cx.set_xlabel("frequency [MHz]",fontsize=fsz)
    ax.set_xlabel("time [$\mu$sec]",fontsize=fsz)
    ax.legend()
    bx.legend()
    cx.legend()
    dx1.legend()
    #dx2.set_ylim([-1,1])

    bx.set_ylabel("Fourier amplitude")
    cx.set_ylabel("group delay [$\mu$sec]")
    dx1.set_ylabel("velocity [km/s]");
    dx2.set_ylabel("inverse decay factor");
    ex.set_ylabel("delay phase angle[deg] ");

    plt.show()
