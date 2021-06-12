import numpy as np
import matplotlib.pyplot as plt

class cp_data:

    def __init__(self):
        fp=open("cp.dat","r")

        fp.readline()
        dat=fp.readline().strip(); dat=dat.split(" ")
        Xa=[float(dat[0]),float(dat[1])]
        print(Xa)

        dat=fp.readline().strip(); dat=dat.split(" ")
        Xb=[float(dat[0]),float(dat[1])]

        fp.readline()
        dat=fp.readline().strip(); dat=dat.split(" ")
        Ndiv=[int(dat[0]),int(dat[1])]
        fp.readline()

        print("Ndiv=",Ndiv)
        cp=[]
        for row in fp:
            cp.append(float(row))
        cp=np.array(cp)
        cp=np.reshape(cp,Ndiv)
        print(np.shape(cp))
        cp=np.transpose(cp)
        
        self.Xa=Xa
        self.Xb=Xb
        self.Ndiv=Ndiv
        self.cp=cp
    def show_cp(self,ax):
        Xa=self.Xa
        Xb=self.Xb
        ext=[Xa[0],Xb[0],Xa[1],Xb[1]]
        im=ax.imshow(self.cp,extent=ext,cmap="jet",origin="lower",interpolation="none")
        return(im);


if __name__=="__main__":

    Vel=cp_data()

    fig=plt.figure()
    ax=fig.add_subplot(111)

    im=Vel.show_cp(ax)
    plt.colorbar(im)
    plt.show()
