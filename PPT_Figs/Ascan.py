# -*- coding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt


class AWV:
    def blank(self,ndat):
        Nt=ndat
        t0=0.0
        amp=np.zeros(ndat)
        time=np.arange(ndat)
        dt=1
    def set_taxis(self,t1,t2,ndat=0):
        if ndat != 0:
            Nt=ndat
        time=np.linspace(t1,t2,Nt)
        dt=time[1]-time[0]
        df=1/(dt*Nt)
        freq=df*np.arange(Nt)

    def load(self,fname):
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
        self.df=df
        self.Nt=len(amp)

    def FFT(self):
        Amp=np.fft.fft(self.amp);
        self.Amp=Amp
    def Butterworth(self,tb,w_6dB,mexp,apply=True):
        tt=self.time-tb;
        Phi=1+(tt/w_6dB)**mexp
        self.Phi=1/Phi
        self.amp*=self.Phi
    def plot_ascan(self,ax):
        ax.plot(self.time,self.amp)
        ax.grid(True)
    def plot_FFT(self,ax,FFT=True):
        if FFT:
            self.Amp=np.fft.fft(self.amp)
        ax.semilogy(self.freq,abs(self.Amp))
        ax.set_xlim([self.freq[1],3])
        ax.set_ylim([1.e-03,10])
        ax.grid(True)
    def Wiener(self,eps=0.2,apply=True,FFT=True):
        if FFT:
            AWV.FFT(self)
        Smax=np.max(self.Amp)
        S2=abs(self.Amp)/Smax
        S2*=S2
        Wnr=S2/(S2+eps*eps)
        self.Amp*=Wnr
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


    def plot_gdelay(self,ax,FFT=True):
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

        ax.plot(self.freq,tg,".",markersize=3)
        ax.set_xlim([0,3])
        ax.grid(True)

def file_name(Mineral,num):

    dir_name="../"
    if Mineral=="Qt":
        dir_name=dir_name+"Quartz"
    if Mineral=="Na":
        dir_name=dir_name+"Na_Feldspar"
    if Mineral=="K":
        dir_name=dir_name+"K_Feldspar"
    fname=dir_name+"/scope_"+str(num)+".csv"

    if Mineral=="Ref":
        fname=dir_name+"1MHznew.csv"

    return(fname)


if __name__=="__main__":

    awv=AWV()

    nums=[111,232,13,13]
    fsz=14
    M=["Qt","Na","K","Ref"]
    clrs=["k","b","r","g"]
    axes=[]
    figs=[]
    for k in range(len(M)):
        Mineral=M[k]
        fname=file_name(Mineral,nums[k])
        clr=clrs[k]
        fig=plt.figure(figsize=[10,4])
        print(fname)
        awv.load(fname)
        ax=fig.add_subplot(111)
        ax.plot(awv.time,awv.amp,"-"+clr,linewidth=2)
        ax.set_xlim([5,30])
        ax.grid(True)
        axes.append(ax)
        figs.append(fig)

        ax.tick_params(labelsize=fsz)
        ax.set_xlabel("time [\mu sec]",fontsize=fsz);
        ax.set_ylabel("amp. [m/s]",fontsize=fsz);

        fig.savefig("typical_"+Mineral+".png",bbox_inches="tight")


    """
    num=300
    dir_name="../"
    if Mineral=="Na":
        dir_name="Na_Feldspar"
        tb=12.5; tw_6dB=1; mexp=6

    if Mineral=="K":
        dir_name="K_Feldspar"
        tb=13; tw_6dB=1; mexp=6

    if Mineral=="Qt":
        dir_name="Quartz"
        tb=13; tw_6dB=1; mexp=6

    fname=dir_name+"/scope_"+str(num)+".csv"

    if Mineral=="Ref":
        fname=dir_name+"1MHznew.csv"
        tb=12; tw_6dB=1; mexp=6
        tb=12; tw_6dB=2.0; mexp=6

    fig2=plt.figure()
    bx=fig2.add_subplot(211)
    cx=fig2.add_subplot(212)

    awv.plot_FFT(bx,FFT=True)
    awv.Butterworth(tb,tw_6dB,mexp,apply=True)
    awv.plot_FFT(bx,FFT=True)
    awv.Wiener(eps=0.2,FFT=False)
    awv.plot_gdelay(cx,FFT=False)

    bx.grid(True)
    cx.grid(True)

    bx.semilogy(awv.freq,abs(awv.Amp))
    awv.plot_ascan(ax)
    """

    plt.show()
