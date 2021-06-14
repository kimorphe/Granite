import matplotlib.pyplot as plt
import numpy as np


class Mmap:
    def load(self,fname):
        fp=open(fname,"r")
        fp.readline()
        npk=int(fp.readline())
        fp.readline()
        dat=fp.readline().strip().split(",")
        nx=int(dat[0])
        ny=int(dat[1])
        print(nx,ny)

        fp.readline()
        M=[];
        for row in fp:
            M.append(int(row))

        M=np.array(M)
        M=np.reshape(M,[nx,ny])
        self.nx=nx;
        self.ny=ny;
        self.M=M
        fp.close()
    def show(self,ax,nclr=7):
        ax.imshow(np.mod(self.M,nclr),aspect=1.0,cmap="hot",vmin=0,vmax=nclr-1)


if __name__=="__main__":
    mp=Mmap()
    mp.load("grain.out")

    fig=plt.figure()
    ax=fig.add_subplot(111)
    mp.show(ax)

    fig.savefig("grain.png",bbox_inches="tight")
    plt.show()


