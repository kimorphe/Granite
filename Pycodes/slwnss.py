# -*- coding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import Ascan as Asc

class Bwv:
    def __init__(self,Nx,Nt):
        self.B=np.zeros([Nx,Nt])
        Amp=np.zeros([Nx,Nt])
        self.Amp=Amp.astype(complex)
        self.tg=np.zeros([Nx,Nt])
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
        return(indx)


class DATA_DIR:
    def __init__(self,Mineral):
        self.Mineral=Mineral
        if Mineral=="Na":
            dir_name="Na_Feldspar"
            tb=12.5; tw_6dB=1; mexp=6
            nfile=548;
            bT=2.03 # Na
        elif Mineral=="K":
            dir_name="K_Feldspar"
            tb=13; tw_6dB=1; mexp=6
            nfile=824;
            bT=1.89 # K
        elif Mineral=="Qt":
            dir_name="Quartz"
            tb=12.5; tw_6dB=1; mexp=6
            nfile=589;
            bT=2.08 #Qt
        else:
            print("no data for ",Mineral)
            exit();
        self.dir_name=dir_name
        self.tb=tb
        self.tw_6dB=tw_6dB
        self.mexp=mexp
        self.nfile=nfile
        self.bT=bT
    def get_dir_name(self):
        return(self.dir_name)
    def get_win_params(self):
        return(self.tb, self.tw_6dB, self.mexp)
    def get_nfile(self):
        return(self.nfile)
    def get_bT(self):
        return(self.bT)
    def show(self):
        print("data directory= ",self.dir_name)
        print("win params(tb, tw, mexp)= ",self.tb,self.tw_6dB,self.mexp)
        print("number of files= ",self.nfile)
        print("efficiency factor= ",self.bT)

if __name__=="__main__":

    Mineral="Na"
    Mineral="K"
    Mineral="Qt"
    ddir=DATA_DIR(Mineral)
    ddir.show()
    dir_name=ddir.get_dir_name()
    tg,tw_6dB,mexp=ddir.get_win_params()
    nfile=ddir.get_nfile()
    bT=ddir.get_bT()

    nums=np.arange(0,nfile,1) # File numbers
    nums=nums.astype(int)

    tg_ref=12.0;    # Window parameter (reference)
    tw_ref=0.8;    # Window parameter (reference)

    # Setup Reference Signal
    aref=Asc.AWV()  # A-scan wave container (reference) 
    aref.load("1MHznew.csv") #load reference waveform
    aref.Butterworth(tg_ref,tw_ref,mexp,apply=True) # Apply window
    aref.gdelay(FFT=True) # Evaluate  group delay


    # Setup B-scan wave data
    awv=Asc.AWV() # A-scan wave container (signal) 
    nfile=len(nums) # Total number of data files (waveforms)

    isum=0  # init. counter
    for k in nums: 
        fname=dir_name+"/scope_"+str(k)+".csv"  # data file name
        awv.load(fname) # load A-scan data
        awv.Butterworth(tg,tw_6dB,mexp,apply=True) # windowing
        #awv.Wiener(eps=0.01,FFT=True) # apply Wiener filter
        awv.gdelay(FFT=True)   # evaluate group delay
        if isum==0:
            bwv=Bwv(nfile,awv.Nt)   # initialize B-scan container
            bwv.set_time(awv.time[0],awv.time[-1]) # set time & frequency axes
        # Copy from awv --> bwv
        bwv.B[isum,:]+=awv.amp[:]   # time signal
        bwv.Amp[isum,:]+=awv.Amp[:] # frequency spectrum
        bwv.tg[isum,:]+=awv.tg[:]   # group delay
        isum+=1 # increment counter
    fig1=plt.figure()
    ax=fig1.add_subplot(111)
    bwv.show(ax)
    #-----------------------------------------


    for k in range(bwv.Nx):
        bwv.tg[k,:]-=aref.tg;   # difference from the reference 
        bwv.Amp[k,:]*=(bT/aref.Amp) # deconvolution

    Hist_t=[]
    Hist_A=[]
    for k in range(bwv.Nt):
        hist,bins=np.histogram(bwv.tg[:,k],range=(0,2),bins=60)
        Hist_t.append(np.array(hist))

        hist,bins=np.histogram(abs(bwv.Amp[:,k]),range=(0,1),bins=60)
        Hist_A.append(np.array(hist))

    fig2=plt.figure()
    bx=fig2.add_subplot(212)    # Group delay
    cx=fig2.add_subplot(211)    # Amplitude 

    fLim=[0.4,1.6]
    Hist_t=np.transpose(np.array(Hist_t))
    ext=[awv.freq[0],awv.freq[-1],0,2]
    bx.imshow(Hist_t,extent=ext,origin="lower",cmap="jet",interpolation="bicubic")
    bx.set_aspect("auto")
    bx.set_xlim(fLim)
    bx.set_ylim([0,2])

    Hist_A=np.transpose(np.array(Hist_A))
    extA=[awv.freq[0],awv.freq[-1],0.0,1]
    cx.imshow(Hist_A,extent=extA,origin="lower",cmap="jet",interpolation="bicubic",vmin=0,vmax=100)
    cx.set_aspect("auto")
    cx.set_xlim(fLim)
    cx.set_ylim([0.0,1])


    Amp_mean1=np.average(bwv.Amp,axis=0)
    Amp_mean2=np.average(abs(bwv.Amp),axis=0)
    Amp_std=np.std(abs(bwv.Amp),axis=0)

    tg_mean=np.average(bwv.tg,axis=0)
    tg_std=np.std(bwv.tg,axis=0)

    freq=bwv.freq
    lwd=0.6
    bx.plot(freq,tg_mean,"w")
    bx.plot(freq,tg_mean+tg_std,"w--",linewidth=lwd)
    bx.plot(freq,tg_mean-tg_std,"w--",linewidth=lwd)
    #cx.plot(freq,abs(Amp_mean1),"k")
    cx.plot(freq,Amp_mean2,"w")
    cx.plot(freq,Amp_mean2+Amp_std,"w--",linewidth=lwd)
    cx.plot(freq,Amp_mean2-Amp_std,"w--",linewidth=lwd)
    bx.grid(True)
    cx.grid(True)


    clrs=["b","y","m","g","r","k"]
    nclr=len(clrs)

    fs=[0.5,0.6,0.7,0.8,1.0,1.2,1.4,1.5,1.6]


    fig3=plt.figure()
    sct=fig3.add_subplot(111)
    zb=[]
    yb=[]
    k=0
    for f0 in fs:
        nfl=bwv.get_freq(f0)
        nfh=bwv.get_freq(f0)
        Z=bwv.Amp[:,nfl:nfh+1]
        Y=bwv.tg[:,nfl:nfh+1]
        nshp=np.shape(Z)
        Z=np.reshape(Z,[nshp[0]*nshp[1],1])
        Y=np.reshape(Y,[nshp[0]*nshp[1],1])
        zb.append(np.mean(abs(Z)))
        yb.append(np.mean(Y))
        sct.plot(Y,abs(Z),"."+clrs[k%nclr],alpha=0.5,markersize=3,label=str(f0)+"[MHz]")
        W=np.array( [Y[:,0], abs(Z[:,0])])
        C=np.cov(W)
        C[0,1]/=np.sqrt((C[0,0]*C[1,1]))
        C[1,0]=C[0,1]
        print(np.sqrt(C))
        k+=1

    k=0
    for f0 in fs:
        sct.plot(yb[k],zb[k],"d"+clrs[k%nclr],alpha=1.0,markersize=6,label=str(f0)+"[MHz]",markeredgecolor="k")
        k+=1
    sct.legend()
    sct.grid(True)
    sct.set_xlim([0,2])
    sct.set_ylim([0,1])


    fsz=14
    ax.set_xlabel("time [$\mu$sec]",fontsize=fsz)
    ax.set_ylabel("data No.",fontsize=fsz)

    bx.set_xlabel("frequency [MHz]",fontsize=fsz)
    bx.set_ylabel("group delay [$\mu$sec]",fontsize=fsz)

    #cx.set_xlabel("frequency [MHz]",fontsize=fsz)
    cx.set_ylabel("amplitude",fontsize=fsz)

    sct.set_xlabel("group delay [$\mu$sec]",fontsize=fsz)
    sct.set_ylabel("amplitude",fontsize=fsz)

    ax.tick_params(labelsize=12)
    bx.tick_params(labelsize=12)
    cx.tick_params(labelsize=12)
    sct.tick_params(labelsize=12)


    fig1.savefig("bscan_"+Mineral+".png",bbox_inches="tight")
    fig2.savefig("hist_"+Mineral+".png",bbox_inches="tight")
    fig3.savefig("sct_"+Mineral+".png",bbox_inches="tight")
    plt.show()
