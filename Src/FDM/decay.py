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
    def FFT(self,bx="",show=True):
        self.Amp=np.fft.fft(self.amp,axis=1)
        self.df=1/self.time[-1];
        self.freq=np.array(range(self.Nt))*self.df
        freq=self.freq
        xcod=self.xcod
        ext=[freq[0],freq[-1],xcod[0],xcod[-1]]
        if show:
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
        indx=np.argmax(self.amp,axis=1)
        self.tmax=self.time[indx]
    def min_amp(self):
        self.amin=np.min(self.amp,axis=1)
        indx=np.argmin(self.amp,axis=1)
        self.tmin=self.time[indx]

    def get_Amp(self,freq):
        num=np.argmin(np.abs(freq-self.freq))
        return((self.Amp[:,num]))


if __name__=="__main__":

    bwv=Bmean()

    fname="v1stk.out"

    #fig=plt.figure()
    #bx=fig.add_subplot(111)

    fig5=plt.figure()
    ex=fig5.add_subplot(111)
    freqs=[1.0]
    ex.grid(True)

    DIR="DAT"
    DIR="Gss"
    nums=[1,2,3,4,5,6,7,8,9,10]
    isum=0
    for num in nums:
    #for dir_name in dircs:
        fn=DIR+str(num)+"/"+fname
        bwv.load(fn)
        print(fn)

        bwv.min_amp()
        t1=bwv.tmin[0]; 
        t2=bwv.tmin[-1];
        print("t1=",t1)
        print("t2=",t2)
        #bwv.show(bx)
        #bx.plot(bwv.tmin,bwv.xcod)
        bwv.Win(t1,t2,1.0); 
        plt.show()
        bwv.FFT(bx="",show=False)

        for frq in freqs:
            Amp=bwv.get_Amp(frq)
            ex.plot(bwv.xcod, np.abs(Amp))

            if isum==0:
                Asum=np.zeros(len(Amp))*1j
            Asum+=Amp;
            isum+=1

    Asum/=isum;
    ex.plot(bwv.xcod,np.abs(Asum),"k-",linewidth=3)

    fp=open(DIR+"_decay.out","w")
    for k in range(len(Asum)):
        dat=str(bwv.xcod[k])+", "+str(np.abs(Asum[k]))+"\n"
        fp.write(dat)
    fp.close()
    plt.show()

