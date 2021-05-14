import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar

class cp1D:
    def load(self,fname):
        fp=open(fname,"r")
        freq=[]
        cp=[]
        for row in fp:
            dat=row.strip().split(",")
            freq.append(float(dat[0]))
            cp.append(float(dat[1]))
        self.freq=freq
        self.cp=cp
    def plot(self,ax,fsz=14,clr="w",lbl=""):
        ax.plot(self.freq,self.cp,"-"+clr,linewidth=2.0,label=lbl)
        #ax.grid(True)
        #ax.set_xlabel("phase velocity [km/s]",fontsize=fsz)
        #ax.set_ylabel("count",fontsize=fsz)

class Hst2D:
    def load(self,fname):
        fp=open(fname);

        fp.readline()
        tmp=fp.readline()

        dat=tmp[1:-1];
        dat=dat.strip().split(",")
        f1=float(dat[0])
        f2=float(dat[1])
        nf=int(dat[2]);

        fp.readline()
        tmp=fp.readline()

        dat=tmp[1:-1]
        dat=dat.strip().split(",")
        c1=float(dat[0])
        c2=float(dat[1])
        nc=int(dat[2]);

        print(f1,f2,nf)
        print(c1,c2,nc)

        tmp=fp.readline()
        print(tmp)

        Z=[]
        for j in range(nf):
            for k in range(nc):
                dat=fp.readline()
                dat=dat.strip().split(",")
                Z.append(float(dat[1]))
            fp.readline()
        print(len(Z))
        Z=np.array(Z)
        Z=np.reshape(Z,[nf,nc])
        Z=np.transpose(Z)

        self.Z=Z           
        self.c1=c1
        self.c2=c2
        self.f1=f1
        self.f2=f2

    def show(self,ax,fsz=14):
        ext=[self.f1,self.f2,self.c1,self.c2]
        im=ax.imshow(self.Z,extent=ext,interpolation="bicubic",aspect="auto",cmap="jet",origin="lower");
        ax.set_xlabel("frequency [MHz]",fontsize=fsz)
        ax.set_ylabel("phase velocity [km/s]",fontsize=fsz)
        ax.tick_params(labelsize=fsz)
        return(im)

if __name__=="__main__":

    H_Qt=Hst2D()
    H_Na=Hst2D()
    H_K=Hst2D()


    H_Qt.load("cp_hist2d_Qt.out")
    H_Na.load("cp_hist2d_Na.out")
    H_K.load("cp_hist2d_K.out")

    MNRL=["Quartz","Na-Feldspar","K-Feldspar"]


    fig1=plt.figure()
    fig2=plt.figure()
    fig3=plt.figure()
    ax1=fig1.add_subplot(111)
    ax2=fig2.add_subplot(111)
    ax3=fig3.add_subplot(111)

    ax1_div=make_axes_locatable(ax1);
    ax2_div=make_axes_locatable(ax2);
    ax3_div=make_axes_locatable(ax3);
    cax1=ax1_div.append_axes("right",size="7%",pad="2%")
    cax2=ax2_div.append_axes("right",size="7%",pad="2%")
    cax3=ax3_div.append_axes("right",size="7%",pad="2%")


    fsz=14
    ax1.set_title(MNRL[0],fontsize=fsz)
    ax2.set_title(MNRL[1],fontsize=fsz)
    ax3.set_title(MNRL[2],fontsize=fsz)

    im1=H_Qt.show(ax1)
    im2=H_Na.show(ax2)
    im3=H_K.show(ax3)

    cbar1=colorbar(im1,cax=cax1)
    cbar2=colorbar(im2,cax=cax2)
    cbar3=colorbar(im3,cax=cax3)


    fsz=12
    cbar1.ax.tick_params(labelsize=fsz)
    cbar2.ax.tick_params(labelsize=fsz)
    cbar3.ax.tick_params(labelsize=fsz)


    cpw=cp1D()
    Fig=plt.figure()
    Ax=Fig.add_subplot(111)
    Ax.grid(True)
    
    clrs=["k","b","r"]
    M=["Qt","Na","K"]
    axs=[ax1,ax2,ax3]
    for k in range(len(M)):
        mnrl=M[k]
        ax=axs[k]
        fname="cp_mean_"+mnrl+".out"
        cpw.load(fname)
        cpw.plot(ax)
        cpw.plot(Ax,clr=clrs[k],lbl=mnrl)

    Ax.set_ylim([3,8])
    Ax.set_xlim([1,2])
    Ax.tick_params(labelsize=14)
    Ax.set_xlabel("frequency [MHz]",fontsize=fsz)
    Ax.set_ylabel("phase velocity [km/s]",fontsize=fsz)
    Ax.legend(fontsize=fsz)

    Fig.savefig("cp_mean.png",bbox_inches="tight")
    fig1.savefig("hist2d_Qt.png",bbox_inches="tight")
    fig2.savefig("hist2d_Na.png",bbox_inches="tight")
    fig3.savefig("hist2d_K.png",bbox_inches="tight")
    plt.show()


