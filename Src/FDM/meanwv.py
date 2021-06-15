import numpy as np
import matplotlib.pyplot as plt

class Bmean:
    def load(self,fname):
        fp=open(fname,"r")
        fp.readline()
        dat=fp.readline().strip().split(",")
        Xa=float(dat[0])
        dx=float(dat[1])
        Nx=int(dat[2])

        fp.readline()
        dat=fp.readline().strip().split(",")
        dt=float(dat[0])
        Nt=int(dat[1])
        fp.readline()

        Z=np.zeros(Nx*Nt)
        for k in range(Nx*Nt):
            Z[k]=float(fp.readline())
        Z=np.reshape(Z,[Nx,Nt])

        self.amp=Z
        self.xcod= np.array(range(Nx))*dx+Xa
        self.time= np.array(range(Nt))*dt
        self.Nt=Nt
        self.Nx=Nx

        fp.close()
    def Win(self,t1,t2,sig):
        tb=(t2-t1)/(self.Nx-1)*np.arange(self.Nx)+t1

        for k in range(self.Nx):
            arg=(self.time-tb[k])/sig
            arg=-0.5*arg*arg;
            Wt=np.exp(arg)
            self.amp[k,:]*=Wt;

        
    def show(self,ax):
        xcod=self.xcod
        time=self.time
        Z=self.amp
        ext=[time[0],time[-1],xcod[0],xcod[-1]]
        im=ax.imshow(Z,extent=ext,cmap="jet",origin="lower",aspect="auto",vmin=-0.15,vmax=0.15)
        return(im)
    def FFT(self,bx):
        self.Amp=np.fft.fft(self.amp,axis=1)
        self.df=1/self.time[-1];
        self.freq=np.array(range(self.Nt))*self.df
        freq=self.freq
        xcod=self.xcod
        ext=[freq[0],freq[-1],xcod[0],xcod[-1]]
        im=bx.imshow(np.abs(self.Amp*self.Amp),extent=ext,cmap="jet",origin="lower",aspect="auto",interpolation="bilinear",vmin=0,vmax=500)
        bx.set_xlim([0,3])
        return(im)
    def kfplot(self,ax):
        self.Amp=np.fft.fft(self.amp,axis=1)
        self.AMP=np.fft.ifft(self.Amp,axis=0)

        self.df=1/self.time[-1];
        self.dk=1/self.xcod[-1]/(2.*np.pi);

        self.kx=np.array(range(self.Nx))*self.dk
        self.freq=np.array(range(self.Nt))*self.df
        ext=[self.freq[0],self.freq[-1],self.kx[0],self.kx[-1]]
        ax.imshow(np.abs(self.AMP),cmap="jet",extent=ext,aspect="auto",origin="lower",interpolation="bilinear")
        ax.set_xlim([0,3])
        ax.set_ylim([0,0.15])
        ax.grid(True)

    def max_amp(self):
        self.amax=np.max(self.amp,axis=1)


if __name__=="__main__":

    fname="v1stk.out"
    #fname="sstk.out"
    bwv=Bmean()
    bwv.load(fname)

    fig=plt.figure()
    ax=fig.add_subplot(111)
    bwv.Win(1.6,11,1.5);
    bwv.Win(1.5,10.6,2.0); 
    im=bwv.show(ax)

    fig2=plt.figure()
    bx=fig2.add_subplot(111)
    im=bwv.FFT(bx)
    plt.colorbar(im)

    fig3=plt.figure()
    cx=fig3.add_subplot(111)
    cx.grid(True)
    bwv.max_amp()
    cx.plot(bwv.xcod,bwv.amax)

    fig4=plt.figure()
    dx=fig4.add_subplot(111)
    bwv.kfplot(dx)



    plt.show()

