import numpy as np
import matplotlib.pyplot as plt
import kw
import data_dir as ddr 
import Gfh

if __name__=="__main__":

    aref=kw.AWV()  # reference wave container
    dref=ddr.DATA_DIR("Ref")
    tb_ref,tw_ref,mexp=dref.get_win_params()
    aref.load("1MHznew.csv")
    aref.Butterworth(tb_ref,tw_ref,mexp,apply=True) # Apply Butterworth window
    aref.FFT()  # FFT 

    ht=3.42 # [mm] sample thickness
    gw=Gfh.Gint(5.0)
    wd=2    # [mm] source width (wedge tip width)
    awv=kw.AWV()   # A-scan wave container

    Mineral="K" # Mineral type 
    Mineral="Qt" # Mineral type 
    Mineral="Na" # Mineral type 

    #   ---- Transmitted Waveform -------
    ddir=ddr.DATA_DIR(Mineral)
    dir_name=ddir.get_dir_name()
    tb,tw_6dB,mexp=ddir.get_win_params()
    nfile=ddir.get_nfile()


    fig1=plt.figure()
    ax=fig1.add_subplot(111)

    num=201     # File No.
    fname=dir_name+"/scope_"+str(num)+".csv"
    ax.set_title(fname)
    awv.load(fname) # load A-scan waveform
    awv.Butterworth(tb,tw_6dB,mexp,apply=True) # apply Butterworth window
    awv.FFT()  # perform FFT
    Wnr=awv.Wiener(eps=0.2,FFT=False,Apply=True) # Wiener filtering
    awv.Amp/=aref.Amp   # Deconvolution
    awv.phase_vel(ht,FFT=False) # obtain phase velocity
    Gm=gw.Decay(awv.freq,awv.cp,ht,wd,301) # evaluate geometrical decay
    """
    awv.Amp=awv.Amp/Gm # geometrical decay compensation
    Z1=dref.rho*dref.cp # acoustic impedance (wedge)
    Z2=ddir.rho*awv.cp  # acoustic impedance (sample)
    Tcoef=2*Z1/(Z1+Z2)  # transmission coefficient
    B=1.77
    awv.Amp=awv.Amp/Tcoef/B
    """
    awv.Amp[0]=0.0
    awv.decay_factor(ht,FFT=False,bT=1.0) # evaulate decay factor (imaginary slowness)
    Ft=awv.Et(rho=ddir.rho) # evaluate temporal relaxation function
    phi=-np.unwrap(np.angle(awv.Amp))

    Et=np.real(awv.Gt)
    H=np.exp(1j*np.angle(awv.Amp))
    r=np.abs(awv.Amp)

    Amp=np.ones(awv.Nt)
    nf1=awv.get_fnum(0.5)
    nf2=awv.get_fnum(1.5)
    Amp[0:nf1]=0;
    Amp[nf2:awv.Nt+1]=0;
    Amp=Amp.astype(complex)
    ofst=0
    #Amp=Wnr.astype(complex)

    for k in range(10):
        hin=np.fft.ifft(Amp)
        ax.plot(awv.time,np.real(hin)+ofst,label=str(k))
        #Amp*=awv.Amp
        #Amp*=(awv.Amp/(np.abs(awv.Amp)+1.e-09))
        Amp*=(H)
        ofst+=0.1
    ax.grid(True)
    ax.set_xlim([0,40])
    plt.show()


    fig2=plt.figure()
    bx1=fig2.add_subplot(211)
    bx2=fig2.add_subplot(212)

    Ew0=np.fft.fft(Et)*awv.dt

    Nt=awv.Nt; Nt2=int(Nt/2)
    Et[Nt2:Nt+1]=0
    ax.plot(awv.time,Et,label="zeroed")

    Ew=np.fft.fft(Et)*awv.dt
    Etd=np.fft.ifft(Ew)*awv.df*awv.Nt
    ax.plot(awv.time,Etd,".",label="IFFT")
    ax.legend()

    #bx1.plot(awv.freq,abs(Ew0))
    #bx1.plot(awv.freq,abs(Ew))
    bx1.grid(True)
    bx2.grid(True)
    #bx2.plot(awv.freq,np.real(Ew0),label="Re{E0}")
    #bx2.plot(awv.freq,np.imag(Ew0),label="Im{E0}")
    #bx2.plot(awv.freq,np.real(Ew),label="Re{E}")
    #bx2.plot(awv.freq,np.imag(Ew),label="Im{E}")

    s0=np.sqrt(ddir.rho/Ew0)
    cp0=1/np.real(s0)
    s20=np.imag(s0)

    s=np.sqrt(ddir.rho/Ew)
    cp=1/np.real(s)
    s2=np.imag(s)
    bx1.plot(awv.freq,cp,label="cp")
    bx2.plot(awv.freq,s2,label="s2")
    bx1.plot(awv.freq,cp0,label="cp0")
    bx2.plot(awv.freq,s20,label="s20")
    bx1.legend()
    bx2.legend()
    bx2.set_ylim([0,1])

    ax.grid(True)

    fig3=plt.figure()
    cx=fig3.add_subplot(111)

    uin=kw.AWV()
    f0=1.0; Nt=awv.Nt;Td=awv.time[-1];
    nbrst=4
    uin.burst(nbrst,f0,Nt,Td)
    T0=1/f0; tb=T0*nbrst*0.5; nsig=6;sig=nbrst*T0/nsig
    uin.Amod_Gauss(tb,sig)
    uin.Amp=np.fft.ifft(uin.amp)
    #uin.FFT()
    nt2=int(uin.Nt*0.5)
    uin.Amp[nt2:uin.Nt+1]=0
    indx=np.argwhere(s2<0)
    s2[indx]=0
    #s=np.real(s)
    #s=np.ones(len(s))/5

    #kin=s*2*np.pi*awv.freq
    kin=(1/awv.cp+1j*awv.s2)*2*np.pi*awv.freq
    kin[0]=0
    xcod=np.array([ht,20,30,40])
    ofst=0
    for xx in xcod:
        Uw=np.exp(1j*kin*xx)#*((uin.Amp))*2
        Uw[Nt2:Nt+1]=0.0
        ut=np.fft.fft(Uw)#*uin.Nt
        cx.plot(uin.time,np.real(ut)+ofst)
        ofst+=0.1

    cx.grid(True)
    plt.show()

