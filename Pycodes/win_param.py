import kw
import numpy as np
import matplotlib.pyplot as plt
import data_dir as ddr


def sigmoid(x,x90):
    a=np.log(0.9/(1-0.9))/x90
    return(1/(1+np.exp(-a*x)))

if __name__=="__main__":


    fig=plt.figure()
    ax=fig.add_subplot(211)
    bx=fig.add_subplot(212)

    fig2=plt.figure()
    cx=fig2.add_subplot(111)

    fig3=plt.figure()
    dx=fig3.add_subplot(111)

    ax.grid(True)
    bx.grid(True)
    cx.grid(True)
    dx.grid(True)

    aref=kw.AWV()

    fname="Na_Feldspar/scope_110.csv"
    fname="Quartz/scope_120.csv"
    fname="K_Feldspar/scope_82.csv"
    fname="1MHznew.csv"
    aref.load(fname,stashed=False)


    aref.plot_ascan(ax,"original")
    aref.plot_ascan(cx,"original")

    t0=11; x90=1
    Zmd=sigmoid(aref.time-t0,x90)
    aref.amp*=Zmd
    aref.plot_ascan(ax,name="Truncated")

    tb=12   # assumed mean TOF
    Sig2=6  # 2-sigma width
    Sig=Sig2/2
    aref.Gauss(tb,Sig,apply=True)   # Preliminary windowing
    aref.plot_ascan(ax)
    amph=aref.Hilbert(FFT=True) # Hilber transform
    zamp=abs(aref.amp+1j*amph)  # analytic signal
    St=np.cumsum(abs(zamp))     # cummulative envelope
    St/=np.max(St)              # normalization
    t50=aref.time[ np.argmin(abs(St-0.5))]  # get t10
    t10=aref.time[ np.argmin(abs(St-0.1))]  # get t50
    t90=aref.time[ np.argmin(abs(St-0.9))]  # get t90
    tb=t50; # Window parameter (mean)
    sig=t50-t10 # Window parameter (std)

    aref.plot_FFT(dx,FFT=False)
    aref.plot_ascan(ax,name="pre-processed")
    aref.plot_ascan(cx,name="pre-processed")
    ax.plot(aref.time,zamp,"k--",label="envelop")
    bx.plot(aref.time,St,"g")

    ylim=ax.get_ylim()
    ax.vlines(t10,ylim[0],ylim[1],"m")
    ax.vlines(t50,ylim[0],ylim[1],"m")
    ax.vlines(t90,ylim[0],ylim[1],"m")
    ax.set_ylim(ylim)

    ylim=bx.get_ylim()
    bx.vlines(t10,ylim[0],ylim[1],"m")
    bx.vlines(t50,ylim[0],ylim[1],"m")
    bx.vlines(t90,ylim[0],ylim[1],"m")
    bx.set_ylim(ylim)

    aref.Gauss(tb,sig,apply=True)   # fine windowing
    #aref.plot_ascan(ax,name="windowed")
    aref.plot_ascan(cx,name="windowed")
    aref.plot_FFT(dx,FFT=True)

    ax.legend()
    print("tb,sig=",tb,sig)


    ax.set_xlim([8,18])
    bx.set_xlim([8,18])
    cx.set_xlim([8,18])


    fig4=plt.figure()
    dx1=fig4.add_subplot(211)
    dx2=fig4.add_subplot(212)

    Minerals=["K","Na","Qt"]
    alph=1
    for M in Minerals:
        ddir=ddr.DATA_DIR(M)
        nfile=ddir.get_nfile()
        dir_name=ddir.get_dir_name()
        tb=[]
        sig=[]
        t90=[]
        for k in range(nfile):
            fname=dir_name+"/scope_"+str(k)+".csv"
            print(fname)
            aref.load(fname,stashed=False)
            (T50,dT,T90)=aref.WinParam()
            tb.append(T50)
            sig.append(dT)
            t90.append(T90)

        dx1.plot(tb,label="t50: "+dir_name)
        dx1.plot(t90,label="t90: "+dir_name)
        #dx2.plot(sig,label=dir_name)
        dx2.hist(sig,alpha=alph,bins=30,range=[0.6,2])
        alph*=0.6
    dx1.grid(True)
    dx2.grid(True)
    dx1.legend()
    dx2.legend()

    print(aref.WinParam())
    plt.show()
