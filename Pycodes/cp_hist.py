import numpy as np
import matplotlib.pyplot as plt

class CpHist:
    def load(self,fname):
        fp=open(fname,"r")

        fp.readline()
        dat=fp.readline()
        dat=dat.strip().split(",")

        f1=float(dat[0])
        f2=float(dat[1])
        nf=int(dat[2])


        fp.readline()
        dat=fp.readline()
        dat=dat.strip().split(",")
        cp1=float(dat[0])
        cp2=float(dat[1])
        nbin=int(dat[2])

        fp.readline()
        Z=[]
        for row in fp:
            Z.append(int(row))
        Z=np.array(Z)
        Z=np.reshape(Z,[nf,nbin])

        self.Z=np.transpose(Z)
        self.freq=np.linspace(f1,f2,nf)
        self.cp=np.linspace(cp1,cp2,nbin)
        self.nbin=nbin
        self.nf=nf

        fp.close()
    def show(self,ax):
        ext=[self.freq[0], self.freq[-1], self.cp[0], self.cp[-1]]
        ax.imshow(self.Z,origin="lower",extent=ext,interpolation="bicubic",cmap="jet")
        fsz=14
        ax.set_xlabel("frequency [MHz]",fontsize=fsz)
        ax.set_ylabel("phase velocity [km/s]",fontsize=fsz)
        ax.set_aspect("auto")



if __name__=="__main__":

    HT=CpHist()
    HT.load("cp.hist")

    fig=plt.figure()
    ax=fig.add_subplot(111)
    HT.show(ax)
    plt.show()

