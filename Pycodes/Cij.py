import numpy as np
import matplotlib.pyplot as plt

class STIFF:
    def load(self,name):
        self.name=name
        self.CIJ=np.zeros([6,6])
        if name=="Qt":
            print("Quartz")
            self.load_Qt()
        if name=="Na":
            print("Na feldspar")
            self.load_Na()
        if name=="K":
            print("K feldspar")
            self.load_K()
        if name=="Bt":
            print("Biotite")
            self.load_Bt()
    def load_Qt(self): # Trigoal
        CIJ=self.CIJ
        # 1,1 block [GPa]
        CIJ[0][0]=87.55; CIJ[0][1]= 6.07; CIJ[0][2]= 13.3;
        CIJ[1][0]=0.000; CIJ[1][1]=87.55; CIJ[1][2]= 13.3;
        CIJ[2][0]=0.000; CIJ[2][1]= 0.00; CIJ[2][2]=106.8;

        # 1,2 block
        CIJ[0][3]= 17.25; CIJ[0][4]= 0.0; CIJ[0][5]= 0.0;
        CIJ[1][3]=-17.25; CIJ[1][4]= 0.0; CIJ[1][5]= 0.0;
        CIJ[2][3]=  0.0;  CIJ[2][4]= 0.0; CIJ[2][5]= 0.0;

        # 2,2 block
        CIJ[3][3]=57.19; CIJ[3][4]= 0.00; CIJ[3][5]= 0.00;
        CIJ[4][3]= 0.00; CIJ[4][4]=57.19; CIJ[4][5]=17.25;
        CIJ[5][3]= 0.00; CIJ[5][4]= 0.00; CIJ[5][5]=40.74;

        self.rho=2.66 # [g/cm3]

        for I in range(6):
            for J in range(6):
                if I <= J:
                    continue
                CIJ[I,J]=CIJ[J,I]
    def load_Na(self): #Triclinic
        CIJ=self.CIJ
        # 1,1 block
        CIJ[0][0]=74.90; CIJ[0][1]= 36.3; CIJ[0][2]= 37.6;
        CIJ[1][0]=0.000; CIJ[1][1]=137.5; CIJ[1][2]= 32.6;
        CIJ[2][0]=0.000; CIJ[2][1]= 0.00; CIJ[2][2]=128.9;

        # 1,2 block
        CIJ[0][3]= 0.00; CIJ[0][4]= -3.1; CIJ[0][5]= 0.0;
        CIJ[1][3]= 0.00; CIJ[1][4]=-10.4; CIJ[1][5]= 0.0;
        CIJ[2][3]= 0.00; CIJ[2][4]=-19.1; CIJ[2][5]= 0.0;

        # 2,2 block
        CIJ[3][3]= 17.2; CIJ[3][4]= 0.00; CIJ[3][5]= -1.3;
        CIJ[4][3]= 0.00; CIJ[4][4]=30.30; CIJ[4][5]=  0.0;
        CIJ[5][3]= 0.00; CIJ[5][4]= 0.00; CIJ[5][5]= 31.1;

        self.rho=2.61

        for I in range(6):
            for J in range(6):
                if I <= J:
                    continue
                CIJ[I,J]=CIJ[J,I]

    def load_K(self): #Triclinic
        CIJ=self.CIJ
        # 1,1 block
        CIJ[0][0]=66.40; CIJ[0][1]= 43.8; CIJ[0][2]= 25.9;
        CIJ[1][0]=0.000; CIJ[1][1]=171.0; CIJ[1][2]= 19.2;
        CIJ[2][0]=0.000; CIJ[2][1]= 0.00; CIJ[2][2]=121.5;

        # 1,2 block
        CIJ[0][3]= 0.00; CIJ[0][4]= -3.3; CIJ[0][5]= 0.0;
        CIJ[1][3]= 0.00; CIJ[1][4]=-14.8; CIJ[1][5]= 0.0;
        CIJ[2][3]= 0.00; CIJ[2][4]=-13.1; CIJ[2][5]= 0.0;

        # 2,2 block
        CIJ[3][3]= 14.3; CIJ[3][4]= 0.00; CIJ[3][5]= -1.5;
        CIJ[4][3]= 0.00; CIJ[4][4]=23.80; CIJ[4][5]=  0.0;
        CIJ[5][3]= 0.00; CIJ[5][4]= 0.00; CIJ[5][5]= 36.1;

        self.rho=2.56

        for I in range(6):
            for J in range(6):
                if I <= J:
                    continue
                CIJ[I,J]=CIJ[J,I]

    def load_Bt(self): # transversely isotropic (heagonal)
        CIJ=self.CIJ
        # 1,1 block
        CIJ[0][0]=186.0; CIJ[0][1]= 32.4; CIJ[0][2]= 11.6;
        CIJ[1][0]=0.000; CIJ[1][1]=186.0; CIJ[1][2]= 11.6;
        CIJ[2][0]=0.000; CIJ[2][1]= 0.00; CIJ[2][2]= 54.0;

        # 1,2 block
        CIJ[0][3]= 0.00; CIJ[0][4]=  0.0; CIJ[0][5]= 0.0;
        CIJ[1][3]= 0.00; CIJ[1][4]=  0.0; CIJ[1][5]= 0.0;
        CIJ[2][3]= 0.00; CIJ[2][4]=  0.0; CIJ[2][5]= 0.0;

        # 2,2 block
        CIJ[3][3]= 5.80; CIJ[3][4]= 0.00; CIJ[3][5]= 0.0;
        CIJ[4][3]= 0.00; CIJ[4][4]= 5.80; CIJ[4][5]= 0.0;
        CIJ[5][3]= 0.00; CIJ[5][4]= 0.00; CIJ[5][5]=76.8;

        self.rho=3.05

        for I in range(6):
            for J in range(6):
                if I <= J:
                    continue
                CIJ[I,J]=CIJ[J,I]

    def ij2I(self,ii,jj):
        i=ii+1
        j=jj+1
        if ii==jj:
            return(i-1)
        else:
            return(9-i-j-1)

    def Gamma(self,xi):
        Gmm=np.zeros([3,3])
        for i in range(3):
            for j in range(3):
                I=self.ij2I(i,j)
                for k in range(3):
                    for l in range(3):
                        J=self.ij2I(k,l)
                        Gmm[i,j]+=(self.CIJ[I,J]*xi[k]*xi[l])
        self.Gmm=Gmm

    def set_Th(self,th1,th2,nth):
        self.Th=np.radians(np.linspace(th1,th2,nth))
        self.nth=nth;
    def set_Phi(self,phi1,phi2,nphi):
        self.Phi=np.radians(np.linspace(phi1,phi2,nphi))
        self.nphi=nphi
    def phase_vels(self):
        isum=0
        self.Pvel=np.zeros(self.nphi*self.nth)
        self.wgt=np.zeros(self.nphi*self.nth)
        for k in range(self.nphi):
            phi=self.Phi[k]
            cosp=np.cos(phi)
            sinp=np.sin(phi)
            for l in range(self.nth):
                th=self.Th[l]
                cost=np.cos(th)
                sint=np.sin(th)
                xi=[sint*cosp, sint*sinp, cost]
                self.Gamma(xi)
                
                W,V=np.linalg.eig(self.Gmm)
                c=np.sqrt(np.abs(W/self.rho))
                self.Pvel[isum]=np.max(c)
                self.wgt[isum]=sint;
                isum+=1


    def plot_sec(self,ax):

        Pvel=np.reshape(self.Pvel,[self.nphi,self.nth])
        for k in range(self.nphi):
            ax.plot(np.degrees(self.Th),Pvel[k,:])
        ax.grid(True)
        fsz=14
        ax.set_xlabel(r"polar angle $\theta$[deg]",fontsize=fsz)
        ax.set_ylabel("phase velocity [km/s]",fontsize=fsz)
        ax.set_title("Quasi P-wave",fontsize=fsz)

    def hist(self,ax,clr="k"):
        nbin=50
        hist,bins=np.histogram(self.Pvel,weights=self.wgt,bins=nbin,range=(3,9))
        hist,bins=np.histogram(self.Pvel,weights=self.wgt,bins=nbin,range=(3,8))
        Bin=0.5*(bins[0:-1]+bins[1:])
        ax.plot(Bin,hist,"-"+clr,linewidth=2)
        ax.grid(True)

        fsz=14
        ax.set_xlabel("phase velocity [km/s]",fontsize=fsz)
        ax.set_ylabel("frequency",fontsize=fsz)
        ax.tick_params(labelsize=fsz)
        print("Mean=",np.average(self.Pvel,weights=self.wgt));
        print("Min=",np.min(self.Pvel));
        print("Max=",np.max(self.Pvel));
if __name__=="__main__":

    fig1=plt.figure()
    ax=fig1.add_subplot(111)

    C=STIFF(); 
    C.set_Th(0,180,181)
    C.set_Phi(0,360,361)


    M=["Qt","Na","K"]
    clrs=["k","b","r"]
    #M=["Qt","Na","K","Bt"]
    #clrs=["k","b","r","m"]

    for k in range(len(M)):
        mnrl=M[k]
        C.load(mnrl)
        print(C.CIJ)
        C.phase_vels()
        C.hist(ax,clr=clrs[k])

    fig1.savefig("hist_cp_mono.png",bbox_inches="tight")
    plt.show()

    fig2=plt.figure()
    bx=fig2.add_subplot(111)
    #C.plot_sec(bx)
