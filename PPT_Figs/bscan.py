# -*- coding: UTF-8 -*-

import numpy as np
import matplotlib.pyplot as plt
import Ascan as Asc
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar

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
        im=ax.imshow(self.B,extent=ext,cmap="jet",origin="lower",interpolation="none",vmin=-0.04,vmax=0.04)
        ax.set_aspect("auto")
        return(im)
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
            dir_name="../Na_Feldspar"
            tb=12.5; tw_6dB=1; mexp=6
            nfile=548;
            bT=2.03 # Na
        elif Mineral=="K":
            dir_name="../K_Feldspar"
            tb=13; tw_6dB=1; mexp=6
            nfile=824;
            bT=1.89 # K
        elif Mineral=="Qt":
            dir_name="../Quartz"
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

    M=["Qt","Na","K"]
    names=["Quartz","Na Feldspar","K Feldspar"]
    ID=2
    Mineral=M[ID]
    name=names[ID]

    ddir=DATA_DIR(Mineral)
    ddir.show()
    dir_name=ddir.get_dir_name()
    nfile=ddir.get_nfile()

    nums=np.arange(0,nfile,1) # File numbers
    nums=nums.astype(int)

    # Setup B-scan wave data
    awv=Asc.AWV() # A-scan wave container (signal) 
    nfile=len(nums) # Total number of data files (waveforms)

    isum=0  # init. counter
    for k in nums: 
        fname=dir_name+"/scope_"+str(k)+".csv"  # data file name
        awv.load(fname) # load A-scan data
        if isum==0:
            bwv=Bwv(nfile,awv.Nt)   # initialize B-scan container
            bwv.set_time(awv.time[0],awv.time[-1]) # set time & frequency axes
        # Copy from awv --> bwv
        bwv.B[isum,:]+=awv.amp[:]   # time signal
        isum+=1 # increment counter
    fig1=plt.figure()
    ax=fig1.add_subplot(111)
    ax_div=make_axes_locatable(ax)
    cax=ax_div.append_axes("right",size="7%",pad="2%")
    im=bwv.show(ax)
    cbar=colorbar(im,cax=cax)
    cbar.ax.tick_params(labelsize=12)
    #-----------------------------------------


    fsz=14
    ax.set_xlabel("time [$\mu$sec]",fontsize=fsz)
    ax.set_ylabel("data No.",fontsize=fsz)
    ax.set_xlim([5,35])

    ax.tick_params(labelsize=fsz)
    ax.set_title(name,fontsize=fsz)
    fig1.savefig("bscan_"+Mineral+".png",bbox_inches="tight")
    plt.show()
