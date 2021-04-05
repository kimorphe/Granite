# -*- coding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import Ascan as Asc

class Bwv:
    def __init__(self,Nx,Nt):
        self.B=np.zeros([Nx,Nt])
        self.Nx=Nx
        self.Nt=Nt

        self.xcod=np.arange(Nx)
        self.time=np.arange(Nt)
        self.t1=0
        self.t2=Nt-1
    def set_time(self,t1,t2):
        self.t1=t1
        self.t2=t2
        self.dt=(t2-t1)/(self.Nt-1)
        self.time=np.arange(self.Nt)*self.dt+self.t1
        self.df=1/(t2-t1);
        self.freq=np.arange(self.Nt)*self.df
    def show(self,ax):
        ext=[self.t1,self.t2,self.xcod[0],self.xcod[self.Nx-1]]
        ax.imshow(self.B,extent=ext,cmap="jet",origin="lower",interpolation="bicubic")
        ax.set_aspect("auto")
    def FFT(self):
        self.Amp=np.fft.fft(self.B,axis=1)
    def show_FFT(self,ax):
        Bwv.FFT(self)
        ext=[self.freq[0],self.freq[-1],self.xcod[0],self.xcod[self.Nx-1]]
        ax.imshow(np.abs(self.Amp),extent=ext,cmap="jet",origin="lower",interpolation="none")
        ax.set_xlim([0,5])
    def gdelay(self):
        Phi=np.angle(self.Amp)
        Phi=np.unwrap(Phi,axis=1)
        dPhi=np.diff(Phi,axis=1)
        tg=np.zeros([self.Nx,self.Nt])
        tg[:,0:self.Nt-1]+=dPhi[:,:]
        tg[:,1:self.Nt]+=dPhi[:,:]
        tg/=2;
        tg[:,0]*=2
        tg[:,-1]*=2
        dw=self.df*2*np.pi
        tg/=dw
        self.tg=-tg
        print(self.tg[0:200])
    def show_gdelay(self,ax):
        Bwv.gdelay(self)
        ext=[self.freq[0],self.freq[-1],self.xcod[0],self.xcod[self.Nx-1]]
        im=ax.imshow(self.tg,extent=ext,cmap="jet",origin="lower",interpolation="none",vmin=0,vmax=20)
        #im=ax.imshow(-self.tg,extent=ext,cmap="jet",origin="lower",interpolation="none")
        ax.set_xlim([0,5])
        return(im)
    def get_freq(self,fr):
        indx=np.argmin(abs(fr-self.freq))
        print(self.freq[indx])
        return(indx)



if __name__=="__main__":

    dir_name="K_Feldspar"
    tg=12.5; tw_6dB=1.0; mexp=6
    tg_ref=12.0;

    aref=Asc.AWV()
    aref.load("1MHznew.csv")
    aref.Butterworth(tg_ref,tw_6dB,mexp,apply=True)
    aref.gdelay()


    awv=Asc.AWV()

    fig1=plt.figure()
    ax=fig1.add_subplot(111)

    nums=np.arange(0,800,1)
    nums=nums.astype(int)
    nfile=len(nums)
    print(nfile)

    isum=0
    for k in nums: 
        fname=dir_name+"/scope_"+str(k)+".csv"
        awv.load(fname)
        awv.Butterworth(tg,tw_6dB,mexp,apply=True)
        #awv.plot_ascan(ax)
        if isum==0:
            bwv=Bwv(nfile,awv.Nt)
            bwv.set_time(awv.time[0],awv.time[-1])
        bwv.B[isum,:]+=awv.amp[:]
        isum+=1
    bwv.show(ax)
    print(np.size(bwv.B))
    print(np.shape(bwv.B))

    #ax.set_xlim([0,30])

    fig2=plt.figure()
    bx=fig2.add_subplot(211)
    cx=fig2.add_subplot(212)

    bwv.FFT()
    bwv.show_FFT(bx)
    im=bwv.show_gdelay(cx)

    bT=2.08 #Qt
    bT=2.03 # Na
    bT=1.89 # K
    for k in range(bwv.Nx):
        bwv.tg[k,:]-=aref.tg;
        bwv.Amp[k,:]*=(bT/aref.Amp)
    cbar=fig2.colorbar(im,ax=cx,orientation="vertical")
    bx.set_aspect("auto")
    cx.set_aspect("auto")


    Amp_mean=np.average(bwv.Amp,axis=0)
    tg_mean=np.average(bwv.tg,axis=0)
    fig4=plt.figure()
    ax_A=fig4.add_subplot(211)
    ax_t=fig4.add_subplot(212)
    ax_A.plot(bwv.freq,abs(Amp_mean))
    ax_t.plot(bwv.freq,tg_mean)
    ax_A.set_xlim([0,3])
    ax_t.set_xlim([0,3])
    ax_A.grid(True)
    ax_t.grid(True)
    



    bx.grid(True)
    cx.grid(True)

    nfl=bwv.get_freq(0.5)
    nfh=bwv.get_freq(0.5)
    print(nfl,nfh)

    fig3=plt.figure()

    Z=bwv.Amp[:,nfl:nfh+1]
    Y=bwv.tg[:,nfl:nfh+1]
    nshp=np.shape(Z)
    print(nshp)
    Z=np.reshape(Z,[nshp[0]*nshp[1],1])
    Y=np.reshape(Y,[nshp[0]*nshp[1],1])
    sct=fig3.add_subplot(111)
    sct.plot(Y,abs(Z),".",alpha=0.5)
    sct.grid(True)

    plt.show()
