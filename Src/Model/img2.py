from PIL import Image
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable
from mpl_toolkits.axes_grid1.colorbar import colorbar

class IMG:
    def __init__(self,fname):
        IM=np.array(Image.open(fname))
        print(type(IM))
        print(IM.dtype)
        print(IM.shape)
        self.IM=IM;
        self.N=IM.shape
    def draw_rect(self,ax,Xa,Xb,clr="k"):
        xs=[Xa[0],Xb[0],Xb[0],Xa[0],Xa[0]]
        ys=[Xa[1],Xa[1],Xb[1],Xb[1],Xa[1]]
        ax.plot(xs,ys,clr)
    def trim(self,Xa,Xb):
        self.IM=self.IM[Xa[0]:Xb[0],Xa[1]:Xb[1],:]
        self.N=self.IM.shape
        print(self.N)
    def show(self,ax,typ=""):
        self.R=self.IM[:,:,0].astype(float)
        self.G=self.IM[:,:,1].astype(float)
        self.B=self.IM[:,:,2].astype(float)
        self.I=(self.R+self.G+self.B)/3.
        R=self.R
        G=self.G
        B=self.B
        I=self.I
        im_cmp="hot"
        if typ=="":
            ax_im=ax.imshow(self.IM)
            ax.set_title("Original")
        if typ=="R":
            ax_im=ax.imshow(self.R,cmap=im_cmp)
            ax.set_title("R")
        if typ=="G":
            ax_im=ax.imshow(self.G,cmap=im_cmp)
            ax.set_title("G")
        if typ=="B":
            ax_im=ax.imshow(self.B,cmap=im_cmp)
            ax.set_title("B")
        if typ=="mean":
            ax_im=ax.imshow(I,cmap=im_cmp)
            ax.set_title("Mean (brightness)")
        if typ=="dR":
            ax_im=ax.imshow(R-I,cmap=im_cmp,vmin=0)
            ax.set_title("dR")
        if typ=="dG":
            ax_im=ax.imshow(G-I,cmap=im_cmp,vmin=0)
            ax.set_title("dG")
        if typ=="dB":
            ax_im=ax.imshow(B-I,cmap=im_cmp,vmin=0)
            ax.set_title("dB")
        if typ=="Bt":
            ax_im=ax.imshow(-I,cmap=im_cmp,vmin=-50,vmax=-0)
            ax.set_title("Biotite")
        if typ=="K":
            ax_im=ax.imshow(R-G,cmap=im_cmp,vmin=-0,vmax=20)
            ax.set_title("K-Feldspar")
        if typ=="Na":
            ax_im=ax.imshow(B-G,cmap=im_cmp,vmin=-5,vmax=20)
            ax.set_title("Na-Feldspar")
        ax.grid(True)
        return(ax_im)
    def show_hist(self,ax,typ=""):
        R=self.R
        G=self.G
        B=self.B
        R=np.reshape(R,[self.N[0]*self.N[1]])
        G=np.reshape(G,[self.N[0]*self.N[1]])
        B=np.reshape(B,[self.N[0]*self.N[1]])
        I=(R+G+B)/3.0
        nbin=50
        if typ=="":
            ax.hist(R,range=[0,255],bins=nbin,color="r")
            ax.hist(G,range=[0,255],bins=nbin,alpha=0.6,color="g")
            ax.hist(B,range=[0,255],bins=nbin,alpha=0.3,color="b")
        if typ=="R":
            ax.hist(R,range=[0,255],bins=nbin,color="r")
            ax.set_title("R")
        if typ=="G":
            ax.hist(G,range=[0,255],bins=nbin,color="g")
            ax.set_title("G")
        if typ=="B":
            ax.hist(B,range=[0,255],bins=nbin,color="b")
            ax.set_title("B")
        if typ=="dR":
            ax.hist(R-I,range=[-30,30],bins=nbin,color="r")
            ax.set_title("dR")
        if typ=="dG":
            ax.hist(G-I,range=[-30,30],bins=nbin,color="g")
            ax.set_title("dG")
        if typ=="dB":
            ax.hist(B-I,range=[-30,30],bins=nbin,color="b")
            ax.set_title("dB")
        if typ=="mean":
            ax.hist(I,range=[0,225],bins=nbin,color="k")
        if typ=="Na":
            ax.hist(B-R,range=[-30,30],bins=nbin,color="G")
        if typ=="K":
            ax.hist(R-G,range=[-25,25],bins=nbin,color="R")
        ax.grid(True)
    def plot_scatter(self,ax):
        R=np.reshape(self.R,[self.N[0]*self.N[1]])
        G=np.reshape(self.G,[self.N[0]*self.N[1]])
        B=np.reshape(self.B,[self.N[0]*self.N[1]])
        I=(R+G+B)/3.0
        #ax.plot(I,R-B,"o",markersize=1,alpha=0.3)
        ax.grid(True)
        ax.set_xlabel("mean (brightness)")
        ax.set_ylabel("reddishness R-B")

        Y2=30
        Y1=-30
        W=np.zeros((255,Y2-Y1))
        for k in range(len(I)):
            ix=int(I[k])
            #iy=int(R[k]-B[k]-Y1)
            iy=int(R[k]-I[k]-Y1)
            #iy=int(I[k]-B[k]-Y1)
            if iy>=Y2-Y1:
                continue
            if iy<0:
                continue
            W[ix,iy]+=1
        W=np.transpose(W)
        ext=[0,225,Y1,Y2]
        ax.imshow(W,cmap="jet",origin="lower",extent=ext,interpolation="bicubic",aspect="auto")

    def min_map(self):
        n1=self.N[0]
        n2=self.N[1]
        y1=40
        y2=125
        y2=120
        M=np.zeros((n1,n2))
        for k in range(n1):
            for l in range(n2):
                r=self.R[k,l]
                g=self.G[k,l]
                b=self.B[k,l]
                y=(r+g+b)/3.

                if y < y1: # Biotite
                    M[k,l]=0 
                elif y < y2: # Quartz
                    M[k,l]=1
                elif r-g>0: # K-feldspar
                    M[k,l]=2 
                else: # Na-feldspar
                    M[k,l]=3
        self.M=M
        return(M)
    def write_M(self,fname):
        fp=open(fname,"w")
        n1=self.N[0]
        n2=self.N[1]
        fp.write("# Nx, Ny\n");
        fp.write(str(n1)+", "+str(n2)+"\n")
        fp.write("# 0:Bt, 1: Qt, 2: K, 3: Na\n")
        for k in range(n1):
            for l in range(n2):
                dat=str(int(self.M[k,l]))+"\n"
                fp.write(dat)
        fp.close();
    def show_min_map(self,ax):
        M=self.M
        R4=np.zeros(np.shape(M))
        G4=np.zeros(np.shape(M))
        B4=np.zeros(np.shape(M))

        # Biotite
        indx=np.where(M==0)
        r=int(np.mean(self.R[indx]))
        g=int(np.mean(self.G[indx]))
        b=int(np.mean(self.B[indx]))
        y=(r+g+b)/3.
        self.Bt_rgb=[r,g,b,y]
        R4[indx]=r; G4[indx]=g; B4[indx]=b;

        # Quartz 
        indx=np.where(M==1)
        r=int(np.mean(self.R[indx]))
        g=int(np.mean(self.G[indx]))
        b=int(np.mean(self.B[indx]))
        y=(r+g+b)/3.
        self.Qt_rgb=[r,g,b,y]
        R4[indx]=r; G4[indx]=g; B4[indx]=b;

        # K-Feldspar 
        indx=np.where(M==2)
        r=int(np.mean(self.R[indx]))
        g=int(np.mean(self.G[indx]))
        b=int(np.mean(self.B[indx]))
        y=(r+g+b)/3.
        self.K_rgb=[r,g,b,y]
        R4[indx]=r; G4[indx]=g; B4[indx]=b;

        # Na-Feldspar 
        indx=np.where(M==3)
        r=int(np.mean(self.R[indx]))
        g=int(np.mean(self.G[indx]))
        b=int(np.mean(self.B[indx]))
        y=(r+g+b)/3.
        self.Na_rgb=[r,g,b,y]
        R4[indx]=r; G4[indx]=g; B4[indx]=b;
        
        M4=self.IM.copy()
        M4[:,:,0]=R4
        M4[:,:,1]=G4
        M4[:,:,2]=B4
        M4img=Image.fromarray(M4)
        ax.imshow(M4img)
        ax.grid(True)

        self.M4img=M4img

if __name__=="__main__":

    #-----------RAW IMAGE---------------
    fig=plt.figure()
    ax=fig.add_subplot(111)
    fname="core_top_face.png"
    img=IMG(fname)
    img.show(ax)

    #-----------TRIMMED IMAGE-----------
    Xa=[200,200]
    Xb=[1000,1000]
    Xa=[300,300]
    Xb=[900,900]
    img.draw_rect(ax,Xa,Xb)
    img.trim(Xa,Xb)

    fig1=plt.figure()
    ax1=fig1.add_subplot(121)
    ax2=fig1.add_subplot(122)
    img.show(ax1)
    img.show(ax2,typ="mean")

    #-----------RGB components-----------
    fig2=plt.figure()
    bx1=fig2.add_subplot(321)
    bx2=fig2.add_subplot(323)
    bx3=fig2.add_subplot(325)
    bx4=fig2.add_subplot(322)
    bx5=fig2.add_subplot(324)
    bx6=fig2.add_subplot(326)

    im1=img.show(bx1,typ="dR")
    im2=img.show(bx2,typ="dG")
    im3=img.show(bx3,typ="dB")
    bx1_divider=make_axes_locatable(bx1)
    bx2_divider=make_axes_locatable(bx2)
    bx3_divider=make_axes_locatable(bx3)
    cax1=bx1_divider.append_axes("right",size="7%",pad="2%")
    cax2=bx2_divider.append_axes("right",size="7%",pad="2%")
    cax3=bx3_divider.append_axes("right",size="7%",pad="2%")
    cbar1=colorbar(im1,cax=cax1)
    cbar2=colorbar(im2,cax=cax2)
    cbar3=colorbar(im3,cax=cax3)

    img.show_hist(bx4,typ="dR")
    img.show_hist(bx5,typ="dG")
    img.show_hist(bx6,typ="dB")

    #-----------Scatter Plot -----------
    fig3=plt.figure()
    cx=fig3.add_subplot(111)
    img.plot_scatter(cx)

    #-----------Mineral Map -----------
    M=img.min_map()
    fig4=plt.figure()
    ex1=fig4.add_subplot(121)
    ex2=fig4.add_subplot(122)
    img.show(ex1)
    img.show_min_map(ex2)

    rgb=img.Bt_rgb
    cx.plot(rgb[3],rgb[0]-rgb[2],"ok",label="Bt")
    rgb=img.Qt_rgb
    cx.plot(rgb[3],rgb[0]-rgb[2],"og",label="Qt")
    rgb=img.Na_rgb
    cx.plot(rgb[3],rgb[0]-rgb[2],"ow",label="Na")
    rgb=img.K_rgb
    cx.plot(rgb[3],rgb[0]-rgb[2],"or",label="K")
    cx.legend()

    img.write_M("minmap.out")

    plt.show()

"""
fp.write(str(len(R))+"\n");
for k in range(len(R)):
    dat=str(R[k])+", "+str(G[k])+", "+str(B[k])+"\n"
    fp.write(dat)
fp.close()
plt.show()
"""
