import numpy as np
import matplotlib.pyplot as plt
import kw 
import data_dir as ddr 
import Gfh
from copy import deepcopy

if __name__=="__main__":

    #   ---- Reference Waveform ------
    aref=kw.AWV()  # reference wave container
    dref=ddr.DATA_DIR("Ref")
    dref.show()
    tb_ref,tw_ref,mexp=dref.get_win_params()
    fnref="1MHznew.csv"
    aref.load(fnref)
    
    (t50,sig,t90)=aref.WinParam()
    aref.Gauss(t50,sig,apply=True)
    #aref.Butterworth(tb_ref,tw_ref,mexp,apply=True) # Apply Butterworth window
    aref.FFT()  # FFT 

    ht=3.42 # [mm] sample thickness
    gw=Gfh.Gint(5.0)
    wd=2    # [mm] source width (wedge tip width)
    awv=kw.AWV()   # A-scan wave container

    Mineral="Qt" # Mineral type 
    Mineral="Na" # Mineral type 
    Mineral="K" # Mineral type 

    #   ---- Transmitted Waveform -------
    ddir=ddr.DATA_DIR(Mineral)
    ddir.show()
    dir_name=ddir.get_dir_name()
    tb,tw_6dB,mexp=ddir.get_win_params()
    nfile=ddir.get_nfile()

    fig2=plt.figure()
    bx1=fig2.add_subplot(211)
    bx2=fig2.add_subplot(212)

    CP=[]
    S2=[]
    AW=[]
    std_c=[]; std_s=[]; 
    ave_c=[]; ave_s=[]; 
    nums=range(nfile)
    fmin=0.5; fmax=1.8
    nf1=aref.get_fnum(fmin)
    nf2=aref.get_fnum(fmax)
    fmin=aref.freq[nf1]
    fmax=aref.freq[nf2]
    for num in nums:
        fname=dir_name+"/scope_"+str(num)+".csv"
        print(fname)
        awv.load(fname) # load A-scan waveform 
        #awv.Butterworth(tb,tw_6dB,mexp,apply=True) # apply Butterworth window
        (t50,sig,t90)=awv.WinParam()
        awv.Gauss(t50,sig,apply=True)
        awv.FFT()  # perform FFT 
        #awv.Wiener(eps=0.2,FFT=False,Apply=True) # Wiener filtering
        awv.Amp/=aref.Amp   # Deconvolution
        awv.phase_vel(ht,FFT=False) # obtain phase velocity
        #if num==0:
        #    Gm=gw.Decay(awv.freq,awv.cp,ht,wd,301) # evaluate geometrical decay 
        #    Gm=np.mean(Gm[nf1:nf2])
        #awv.Amp=awv.Amp/Gm # geometrical decay compensation

        Z1=dref.rho*dref.cp # acoustic impedance (wedge)
        Z2=ddir.rho*awv.cp  # acoustic impedance (sample)
        Tcoef=2*Z1/(Z1+Z2)  # transmission coefficient
        B=1.77
        awv.Amp=awv.Amp/Tcoef/B
        awv.Amp[0]=0.0
        awv.decay_factor(ht,FFT=False,bT=1.0) # evaulate decay factor (imaginary slowness)
        #awv.gdelay(FFT=False)
        #cg=ht/awv.tg


        #bx1.plot(awv.freq[nf1:nf2],awv.cp[nf1:nf2],"o",markersize=2,alpha=0.8)
        #bx2.plot(awv.freq[nf1:nf2],awv.s2[nf1:nf2]/awv.cp[nf1:nf2],"o",markersize=2,alpha=0.8)

        std_c.append(np.std(awv.cp[nf1:nf2]))
        ave_c.append(np.mean(awv.cp[nf1:nf2]))
        std_s.append(np.std(awv.s2[nf1:nf2]*awv.cp[nf1:nf2]))
        ave_s.append(np.mean(awv.s2[nf1:nf2]*awv.cp[nf1:nf2]))

        CP.append(awv.cp[nf1:nf2])
        S2.append(awv.s2[nf1:nf2])
        AW.append(abs(awv.Amp[nf1:nf2]))

    nfile=len(nums)
    npnt=nf2-nf1
    n=nfile*npnt

    Hist_cp=[]
    Hist_s2=[]
    CP=np.array(CP)
    S2=np.array(S2)
    AW=np.array(AW)
    indx=np.where(CP<0)
    jndx=np.where(CP>8)
    #AW=np.ones(np.shape(CP))
    AW[indx]=0; AW[jndx]=0
    cpw=np.average(CP,axis=0,weights=AW)
    s2w=np.average(S2/CP,axis=0,weights=AW)
    for k in range(npnt):
        hist,bins=np.histogram(CP[:,k],range=(3,8),bins=30)
        Hist_cp.append(np.array(hist))
    ext=[fmin,fmax,3,8]
    Hist_cp=np.transpose(np.array(Hist_cp))
    bx1.imshow(Hist_cp,extent=ext,origin="lower",cmap="jet",interpolation="bicubic")
    bx1.set_aspect("auto")
    bx1.plot(awv.freq[nf1:nf2],cpw,"w--")
    I=np.argmax(Hist_cp,axis=0)
    c0=bins[I]+(bins[1]-bins[0])*0.5
    #bx1.plot(awv.freq[nf1:nf2],c0,"k.")

    for k in range(npnt):
        hist,bins=np.histogram(S2[:,k]/CP[:,k],range=(-0.00,0.05),bins=40)
        Hist_s2.append(np.array(hist))
    ext=[fmin,fmax,0.00,0.05]
    Hist_s2=np.transpose(np.array(Hist_s2))
    bx2.imshow(Hist_s2,extent=ext,origin="lower",cmap="jet",interpolation="bicubic")
    bx2.set_aspect("auto")
    bx2.plot(awv.freq[nf1:nf2],s2w,"w--")
    I=np.argmax(Hist_s2,axis=0)
    s20=bins[I]+(bins[1]-bins[0])*0.5
    #bx2.plot(awv.freq[nf1:nf2],s20,"k.")

    fig5=plt.figure()
    ex1=fig5.add_subplot(211)
    ex2=fig5.add_subplot(212)
    Hist_Ew=[]
    Hist_phi=[]
    Ew=1/CP+1j*S2
    Ew=ddir.rho/(Ew*Ew)
    for k in range(npnt):
        hist,bins=np.histogram(abs(Ew[:,k]),range=(0,150),bins=30)
        Hist_Ew.append(np.array(hist))

        hist,bins=np.histogram(np.degrees(np.angle(Ew[:,k])),range=(-100,40),bins=40)
        Hist_phi.append(np.array(hist))

    ext=[fmin,fmax,0,120]
    Hist_Ew=np.transpose(np.array(Hist_Ew))
    ex1.imshow(Hist_Ew,extent=ext,origin="lower",cmap="jet",interpolation="bicubic")
    ex1.set_aspect("auto")
    ex1.grid(True)
    ex1.set_title(dir_name)

    ext=[fmin,fmax,-100,40]
    Hist_phi=np.transpose(np.array(Hist_phi))
    ex2.imshow(Hist_phi,extent=ext,origin="lower",cmap="jet",interpolation="bicubic")
    ex2.set_aspect("auto")
    ex2.grid(True)

    CP=np.reshape(CP,[n,1])
    S2=np.reshape(S2,[n,1])
    bx1.grid(True)
    bx2.grid(True)
    bx1.set_ylim([3,8])
    bx1.set_title(dir_name)

    fig3=plt.figure()
    cx1=fig3.add_subplot(121)
    cx2=fig3.add_subplot(122)
    cx1.hist(CP,bins=60,range=[3,8])
    cx2.hist(S2/CP,bins=60,range=[-0.02,0.05])
    cx1.grid(True)
    cx2.grid(True)
    bx1.set_title(dir_name)

    fig4=plt.figure()
    dx1=fig4.add_subplot(211)
    dx2=fig4.add_subplot(212)
    ave_c=np.array(ave_c)
    std_c=np.array(std_c)
    ave_s=np.array(ave_s)
    std_s=np.array(std_s)

    dx1.hist(std_c/ave_c,bins=40,label="std{cp}",range=(0,0.5))
    dx2.hist(std_s/ave_s,bins=40,label="std{s2}",range=(0,0.5))
    dx1.grid(True)
    dx2.grid(True)
    dx1.legend()
    dx2.legend()

    k=0
    awv_c0=kw.AWV()
    awv_s20=kw.AWV()
    for num in nums:
        fname=dir_name+"/scope_"+str(num)+".csv"
        print(fname)
        awv.load(fname) # load A-scan waveform 
        (t50,sig,t90)=awv.WinParam()
        awv.Gauss(t50,sig,apply=True)
        awv.FFT()  # perform FFT 
        awv.Amp/=aref.Amp   # Deconvolution
        awv.phase_vel(ht,FFT=False) # obtain phase velocity
        awv.Amp[0]=0.0
        awv.decay_factor(ht,FFT=False,bT=1.0) # evaulate decay factor (imaginary slowness)
        r2=awv.cp[nf1:nf2]-c0
        r2=np.sum(r2*r2)

        q2=awv.s2[nf1:nf2]/awv.cp[nf1:nf2]-s20
        q2=np.sum(q2*q2)
        if k==0:
            r2min=r2
            q2min=q2
            imin=num
            jmin=num
            awv_c0=deepcopy(awv) 
            awv_s20=deepcopy(awv) 
        if r2 < r2min:
            r2min=r2
            imin=num
            awv_c0=deepcopy(awv) 
        if q2 < q2min:
            q2min=q2
            jmin=num
            awv_s20=deepcopy(awv) 
        k+=1

    print(imin,r2min)
    print(jmin,q2min)

    bx1.plot(awv_c0.freq[nf1:nf2],awv_c0.cp[nf1:nf2],"w-")
    bx2.plot(awv_s20.freq[nf1:nf2],awv_s20.s2[nf1:nf2]/awv_s20.cp[nf1:nf2],"w-")
    fig2.savefig(Mineral+"_cs.png",bbox_inches="tight")
    fig5.savefig(Mineral+"_Ephi.png",bbox_inches="tight")
    fig3.savefig(Mineral+"_hist.png",bbox_inches="tight")
    plt.show()
