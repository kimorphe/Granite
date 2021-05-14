import numpy as np
import matplotlib.pyplot as plt

class Hst1D:
    def load(self,fname):
        fp=open(fname,"r")
        cp=[]
        ht=[]
        for row in fp:
            dat=row.strip().split(",")
            cp.append(float(dat[0]))
            ht.append(float(dat[1]))
        self.cp=cp
        self.ht=ht
    def plot(self,ax,fsz=14,clr="k"):
        ax.plot(self.cp,self.ht,"-"+clr,linewidth=2.0)
        ax.grid(True)
        ax.set_xlabel("phase velocity [km/s]",fontsize=fsz)
        ax.set_ylabel("count",fontsize=fsz)


if __name__=="__main__":

    fig=plt.figure()
    ax=fig.add_subplot(111)
    M=["Qt","Na","K"]
    #M=["K","Na","Qt"]

    clrs=["k","b","r"]
    k=0
    for mnrl in M:
        fname="cp_hist_"+mnrl+".out"
        H=Hst1D()
        H.load(fname)
        H.plot(ax,clr=clrs[k])
        k+=1

    ax.tick_params(labelsize=14)
    fig.savefig("hist_total.png",bbox_inches="tight")
    plt.show()


