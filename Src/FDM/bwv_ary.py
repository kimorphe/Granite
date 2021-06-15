import numpy as np
import matplotlib.pyplot as plt


class Bwv:
    def load(self,fname):
        fp=open(fname,"r")
        fp.readline()
        idir=int(fp.readline())
        fp.readline()
        npnt=int(fp.readline())

        fp.readline()
        dat=fp.readline().strip().split(",")
        dt=float(dat[0])
        Nt=int(dat[1])


        v1=[]
        v2=[];
        s=[]
        xcod=[]
        ycod=[]
        for k in range(npnt):
            fp.readline()
            dat=fp.readline().strip().split(",")
            xf=float(dat[0]);
            yf=float(dat[1]);
            #print(xf,yf)

            xcod.append(xf)
            ycod.append(yf)
            for l in range(Nt):
                dat=fp.readline().strip().split(",")
                v1.append(float(dat[0]))
                v2.append(float(dat[1]))
                s.append(float(dat[2]))


        v1=np.array(v1)
        v2=np.array(v2)
        s=np.array(s)

        v1=np.reshape(v1,[npnt,Nt])
        v2=np.reshape(v2,[npnt,Nt])
        s=np.reshape(s,[npnt,Nt])

        self.idir=idir
        self.dt=dt
        self.Nt=Nt
        self.npnt=npnt
        self.v1=v1
        self.v2=v2
        self.s=s

        self.time=np.array(range(Nt))*dt
        self.xcod=np.array(xcod)
        self.ycod=np.array(ycod)

        fp.close()
    def show(self,ax):
        time=self.time
        xcod=self.xcod
        ycod=self.ycod
        ext=[time[0],time[-1],xcod[0],xcod[-1]]
        if self.idir==1:
            ext=[time[0],time[-1],ycod[0],ycod[-1]]
        ax.imshow(self.v1,extent=ext,origin="lower",cmap="jet",interpolation="bilinear")
        ax.set_aspect("auto")

if __name__=="__main__":

    fig=plt.figure()
    ax=fig.add_subplot(111)

    bwv0=Bwv()
    bwv0.load("bwv1.out")
    bwv0.show(ax)


    fig2=plt.figure()
    ax2=fig2.add_subplot(111)
    ax2.plot(bwv0.time,bwv0.v1[0,:])
    ax2.plot(bwv0.time,bwv0.v1[20,:])
    ax2.plot(bwv0.time,bwv0.v1[-1,:])
    ax2.grid(True)
    plt.show()

