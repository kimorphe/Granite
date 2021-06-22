import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar

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
        im=ax.imshow(self.cp,extent=ext,cmap="jet",origin="lower",interpolation="none",vmin=4,vmax=7)
        return(im);


if __name__=="__main__":

    Vel=cp_data()

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax_divider=make_axes_locatable(ax)
    cax=ax_divider.append_axes("right",size="7%",pad="2%");

    im=Vel.show_cp(ax)
    cb=colorbar(im,cax=cax,orientation="vertical")
    #plt.colorbar(im)
    ax.tick_params(labelsize=14)
    cb.ax.tick_params(labelsize=12)
    plt.show()

    fig.savefig("cp.png",bbox_inches="tight")
