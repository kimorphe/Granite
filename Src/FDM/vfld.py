import numpy as np
import matplotlib.pyplot as plt

class vfld:

    def load(self,num):
        fname="v"+str(num)+".out"
        fp=open(fname,"r");

        fp.readline()
        dat=fp.readline().strip(); dat=dat.split(",")
        Xa=[float(dat[0]),float(dat[1])]
        Xb=[float(dat[2]),float(dat[3])]
        print(Xa)
        print(Xb)

        fp.readline()
        dat=fp.readline().strip(); dat=dat.split(",")
        Ndiv=[int(dat[0]),int(dat[1])]
        fp.readline()

        print("Ndiv=",Ndiv)
        vx=[]; vy=[];
        for row in fp:
            dat=row.strip().split(",")
            vx.append(float(dat[0]))
            vy.append(float(dat[1]))
        vx=np.array(vx)
        vy=np.array(vy)
        vx=np.reshape(vx,Ndiv)
        vy=np.reshape(vy,Ndiv)
        vx=np.transpose(vx)
        vy=np.transpose(vy)
        
        self.Xa=Xa
        self.Xb=Xb
        self.Ndiv=Ndiv
        self.vx=vx
        self.vy=vy
    def show_v(self,ax):
        Xa=self.Xa
        Xb=self.Xb
        ext=[Xa[0],Xb[0],Xa[1],Xb[1]]
        Z=np.abs(self.vx*self.vx+self.vy*self.vy)
        v1=0.0
        v2=0.1
        ax.imshow(Z,extent=ext,cmap="jet",origin="lower",interpolation="bilinear",vmin=v1,vmax=v2)


if __name__=="__main__":

    fig=plt.figure()
    ax=fig.add_subplot(111)

    vf=vfld()


    nums=[0,1,2,3,4];

    for num in nums:
        vf.load(num)
        vf.show_v(ax)
        fout="v"+str(num)+".png"
        fig.savefig(fout,bbox_inches="tight")
    #plt.show()
