import numpy as np
import matplotlib.pyplot as plt
import kw 
import data_dir as ddr 
import Gfh
from copy import deepcopy
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar

if __name__=="__main__":


    #   FIGURES & AXES 
    fsz=12
    fig1=plt.figure(figsize=[7,8])
    bx1=fig1.add_subplot(211)   # phase velocity (cp)
    bx2=fig1.add_subplot(212)   # decay factor s2/cp (inverse Q)
    bx2.set_xlabel("frequency [MHz]",fontsize=fsz)
    bx1.set_ylabel("phase velocity [km/s]",fontsize=fsz)
    bx2.set_ylabel("decay factor",fontsize=fsz)

    bx1_divider=make_axes_locatable(bx1)
    bx2_divider=make_axes_locatable(bx2)
    cax1=bx1_divider.append_axes("right",size="7%",pad="2%")
    cax2=bx2_divider.append_axes("right",size="7%",pad="2%")

    fig2=plt.figure(figsize=[7,8])
    ex1=fig2.add_subplot(211)   # Young's modulus |E| (magnitude)
    ex2=fig2.add_subplot(212)   # phase angle 
    ex2.set_xlabel("frequency [MHz]",fontsize=fsz)
    ex1.set_ylabel("magnitude [GPa]",fontsize=fsz)
    ex2.set_ylabel("argument [deg]",fontsize=fsz)

    ex1_divider=make_axes_locatable(ex1)
    ex2_divider=make_axes_locatable(ex2)
    cax3=ex1_divider.append_axes("right",size="7%",pad="2%")
    cax4=ex2_divider.append_axes("right",size="7%",pad="2%")

    fig3=plt.figure(figsize=[11,5])
    cx1=fig3.add_subplot(121)   # histogram (cp)
    cx2=fig3.add_subplot(122)   # histogram (s2/cp)
    cx1.set_xlabel("phase velocity [km/s]",fontsize=fsz)
    cx2.set_xlabel("decay factor",fontsize=fsz)
    cx1.set_ylabel("frequency",fontsize=fsz)

    fig4=plt.figure()
    dx1=fig4.add_subplot(211) # coefficient of variation (cp)
    dx2=fig4.add_subplot(212) # coefficient of variation (s2/cp)

    axes=[bx1,bx2,ex1,ex2,cx1,cx2,dx1,dx2]
    
    #  Input Data
    Mineral="Na" # Mineral type 
    Mineral="K" # Mineral type 
    Mineral="Qt" # Mineral type 

    cmin=3; cmax=8
    smin=-0.00; smax=0.5
    Emin=0; Emax=120
    amin=-60; amax=30
    fmin=0.5; fmax=1.6
    #-----------------------------

    #   ---- Reference Waveform ------
    aref=kw.AWV()  # reference wave container
    dref=ddr.DATA_DIR("Ref")
    dref.show()
    tb_ref,tw_ref,mexp=dref.get_win_params()
    fnref="1MHznew.csv"
    aref.load(fnref)
    
    (t50,sig,t90)=aref.WinParam()
    t50=11.8
    sig=0.5
    aref.Gauss(t50,sig,apply=True)
    #aref.Butterworth(tb_ref,tw_ref,mexp,apply=True) # Apply Butterworth window
    aref.FFT()  # FFT 

    ht=3.42 # [mm] sample thickness
    gw=Gfh.Gint(5.0)
    wd=2    # [mm] source width (wedge tip width)
    awv=kw.AWV()   # A-scan wave container

    #   ---- Transmitted Waveform -------
    ddir=ddr.DATA_DIR(Mineral)
    ddir.show()
    dir_name=ddir.get_dir_name()
    tb,tw_6dB,mexp=ddir.get_win_params()
    nfile=ddir.get_nfile()

    CP=[]
    S2=[]
    AW=[]
    std_c=[]; std_s=[]; 
    ave_c=[]; ave_s=[]; 
    nums=range(nfile)
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
        sig=0.5
        t50=12.5
        awv.Gauss(t50,sig,apply=True)
        awv.FFT()  # perform FFT 
        #awv.Wiener(eps=0.2,FFT=False,Apply=True) # Wiener filtering
        awv.Amp/=aref.Amp   # Deconvolution
        awv.phase_vel(ht,FFT=False) # obtain phase velocity
        #if num==0:
        #    Gm=gw.Decay(awv.freq,awv.cp,ht,wd,301) # evaluate geometrical decay 
        #    Gm=np.mean(Gm[nf1:nf2])
        #-------------
        Gm=1/2.5
        awv.Amp=awv.Amp/Gm # geometrical decay compensation
        Z1=dref.rho*dref.cp # acoustic impedance (wedge)
        Z2=ddir.rho*awv.cp  # acoustic impedance (sample)
        Tcoef=2*Z1/(Z1+Z2)  # transmission coefficient
        B=1.77
        B=2.00
        awv.Amp=awv.Amp/Tcoef/B
        #-------------
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
    AW[indx]=0; AW[jndx]=0
    cpw=np.average(CP,axis=0,weights=AW)
    s2w=np.average(S2*CP,axis=0,weights=AW)
    for k in range(npnt):
        hist,bins=np.histogram(CP[:,k],range=(cmin,cmax),bins=30)
        Hist_cp.append(np.array(hist))
    ext=[fmin,fmax,cmin,cmax]
    Hist_cp=np.transpose(np.array(Hist_cp))
    im1=bx1.imshow(Hist_cp,extent=ext,origin="lower",cmap="jet",interpolation="bicubic")
    bx1.set_aspect("auto")
    #bx1.plot(awv.freq[nf1:nf2],cpw,"w--")
    I=np.argmax(Hist_cp,axis=0)
    c0=bins[I]+(bins[1]-bins[0])*0.5
    #bx1.plot(awv.freq[nf1:nf2],c0,"k.")

    for k in range(npnt):
        hist,bins=np.histogram(S2[:,k]*CP[:,k],range=(smin,smax),bins=40)
        Hist_s2.append(np.array(hist))
    ext=[fmin,fmax,smin,smax]
    Hist_s2=np.transpose(np.array(Hist_s2))
    im2=bx2.imshow(Hist_s2,extent=ext,origin="lower",cmap="jet",interpolation="bicubic")
    bx2.set_aspect("auto")
    #bx2.plot(awv.freq[nf1:nf2],s2w,"w--")
    I=np.argmax(Hist_s2,axis=0)
    s20=bins[I]+(bins[1]-bins[0])*0.5
    #bx2.plot(awv.freq[nf1:nf2],s20,"k.")

    Hist_Ew=[]
    Hist_phi=[]
    Ew=1/CP+1j*S2
    Ew=ddir.rho/(Ew*Ew)
    for k in range(npnt):
        hist,bins=np.histogram(abs(Ew[:,k]),range=(Emin,Emax),bins=30)
        Hist_Ew.append(np.array(hist))
    Hist_Ew=np.transpose(np.array(Hist_Ew))
    I=np.argmax(Hist_Ew,axis=0)
    E0=bins[I]+(bins[1]-bins[0])*0.5

    ext=[fmin,fmax,Emin,Emax]
    im3=ex1.imshow(Hist_Ew,extent=ext,origin="lower",cmap="jet",interpolation="bicubic")
    ex1.set_aspect("auto")

    for k in range(npnt):
        hist,bins=np.histogram(np.degrees(np.angle(Ew[:,k])),range=(amin,amax),bins=40)
        Hist_phi.append(np.array(hist))
    Hist_phi=np.transpose(np.array(Hist_phi))
    I=np.argmax(Hist_phi,axis=0)
    phi0=bins[I]+(bins[1]-bins[0])*0.5

    ext=[fmin,fmax,amin,amax]
    im4=ex2.imshow(Hist_phi,extent=ext,origin="lower",cmap="jet",interpolation="bicubic")
    ex2.set_aspect("auto")

    CP=np.reshape(CP,[n,1])
    S2=np.reshape(S2,[n,1])
    bx1.set_ylim([cmin,cmax])

    cx1.hist(CP,bins=60,range=[0,10])
    cx2.hist(S2*CP,bins=60,range=[-0.1,0.5])

    ave_c=np.array(ave_c)
    std_c=np.array(std_c)
    ave_s=np.array(ave_s)
    std_s=np.array(std_s)

    dx1.hist(std_c/ave_c,bins=40,label="std{cp}",range=(0,0.6))
    dx2.hist(std_s/ave_s,bins=40,label="std{s2}",range=(0,0.6))
    dx1.legend()
    dx2.legend()

    k=0
    awv_c0=kw.AWV()
    awv_s20=kw.AWV()
    awv_E0=kw.AWV()
    awv_phi0=kw.AWV()
    r2min=0; q2min=0
    p2min=0; w2min=0
    for num in nums:
        fname=dir_name+"/scope_"+str(num)+".csv"
        print(fname)
        awv.load(fname) # load A-scan waveform 
        (t50,sig,t90)=awv.WinParam()
        sig=0.5
        t50=12.5

        awv.Gauss(t50,sig,apply=True)
        awv.FFT()  # perform FFT 
        awv.Amp/=aref.Amp   # Deconvolution
        awv.phase_vel(ht,FFT=False) # obtain phase velocity
        #-------------
        Gm=1/2.5
        awv.Amp=awv.Amp/Gm # geometrical decay compensation
        Z1=dref.rho*dref.cp # acoustic impedance (wedge)
        Z2=ddir.rho*awv.cp  # acoustic impedance (sample)
        Tcoef=2*Z1/(Z1+Z2)  # transmission coefficient
        B=1.77
        B=2.00
        awv.Amp=awv.Amp/Tcoef/B
        #-------------
        awv.Amp[0]=0.0
        awv.decay_factor(ht,FFT=False,bT=1.0) # evaulate decay factor (imaginary slowness)
        r2=awv.cp[nf1:nf2]-c0
        r2=np.sum(r2*r2)

        q2=awv.s2[nf1:nf2]*awv.cp[nf1:nf2]-s20
        q2=np.sum(q2*q2)

        Ew=1/awv.cp+1j*awv.s2
        Ew=ddir.rho/(Ew*Ew)
        Phi=np.degrees(np.angle(Ew))
        Ew=np.abs(Ew)
        p2=Ew[nf1:nf2]-E0
        p2=np.sum(p2*p2)

        w2=Phi[nf1:nf2]-phi0
        w2=np.sum(w2*w2)

        if r2 < r2min or k==0:
            r2min=r2
            imin=num
            awv_c0=deepcopy(awv) 
        if q2 < q2min or k==0:
            q2min=q2
            jmin=num
            awv_s20=deepcopy(awv) 
        if p2 < p2min or k==0:
            p2min=p2
            kmin=num
            awv_E0=deepcopy(awv) 
        if w2 < w2min or k==0:
            w2min=w2
            lmin=num
            awv_phi0=deepcopy(awv) 
        k+=1

    print("typical cp  -->",imin)
    print("typical s/c -->",jmin)
    print("typical |E| -->",kmin)
    print("typical arg{E} -->",lmin)
    bx1.plot(awv_c0.freq[nf1:nf2],awv_c0.cp[nf1:nf2],"w-")
    bx2.plot(awv_s20.freq[nf1:nf2],awv_s20.s2[nf1:nf2]*awv_s20.cp[nf1:nf2],"w-")

    E0=np.abs(1/awv_E0.cp+1j*awv_E0.s2)
    E0=np.abs(ddir.rho/(E0*E0))
    ex1.plot(awv_E0.freq[nf1:nf2],E0[nf1:nf2],"w-")

    phi0=1/awv_phi0.cp+1j*awv_phi0.s2
    phi0=ddir.rho/(phi0*phi0)
    phi0=np.degrees(np.angle(phi0))
    ex2.plot(awv_phi0.freq[nf1:nf2],phi0[nf1:nf2],"w-")

    for axis in axes:
        axis.tick_params(labelsize=fsz)
        axis.grid(True)
    bx1.set_title(dir_name)
    ex1.set_title(dir_name)

    cbar1=colorbar(im1,cax=cax1) 
    cbar2=colorbar(im2,cax=cax2) 
    cbar3=colorbar(im3,cax=cax3) 
    cbar4=colorbar(im4,cax=cax4) 
    cbar1.ax.tick_params(labelsize=fsz)
    cbar2.ax.tick_params(labelsize=fsz)
    cbar3.ax.tick_params(labelsize=fsz)
    cbar4.ax.tick_params(labelsize=fsz)

    Fsz=16
    xtxt=(fmax-fmin)*0.05+fmin
    ytxt=(cmax-cmin)*0.9+cmin
    bx1.text(xtxt,ytxt,"(a)",color="w",fontsize=Fsz)
    ytxt=(smax-smin)*0.9+smin
    bx2.text(xtxt,ytxt,"(b)",color="w",fontsize=Fsz)
    ytxt=(Emax-Emin)*0.9+Emin
    ex1.text(xtxt,ytxt,"(a)",color="w",fontsize=Fsz)
    ytxt=(amax-amin)*0.9+amin
    ex2.text(xtxt,ytxt,"(b)",color="w",fontsize=Fsz)

    fig1.savefig(Mineral+"_cs.png",bbox_inches="tight")
    fig2.savefig(Mineral+"_Ephi.png",bbox_inches="tight")
    fig3.savefig(Mineral+"_hist.png",bbox_inches="tight")
    plt.show()
