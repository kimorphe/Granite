import numpy as np
import matplotlib.pyplot as plt


fname="bwv.out"
fp=open(fname,"r")

fp.readline()
dat=fp.readline().strip().split(",")
xr1=float(dat[0])
xr2=float(dat[1])
nrx=int(dat[2])


fp.readline()
dat=fp.readline().strip().split(",")
yr1=float(dat[0])
yr2=float(dat[1])
nry=int(dat[2])

fp.readline()
dat=fp.readline().strip().split(",")
Nt=int(dat[0])
dt=float(dat[1])

time=np.arange(Nt)*dt

fp.readline()


fig=plt.figure()
ax=fig.add_subplot(111)

v1=[]; v2=[]; s=[]
#Nx=nry*nrx;
Nx=nry;
for m in range(Nx):
    for k in range(Nt):
        dat=fp.readline().strip().split(",")
        v1.append(float(dat[0]))
        v2.append(float(dat[1]))
        s.append(float(dat[2]))

v1=np.array(v1)
v2=np.array(v2)
s=np.array(s) 
v1=np.reshape(v1,[Nx,Nt]);
v2=np.reshape(v2,[Nx,Nt]);
s=np.reshape(s,[Nx,Nt]);

ext=[time[0],time[-1],0,Nx]
ax.imshow(v1,origin="lower",cmap="jet",extent=ext,interpolation="bilinear")
ax.set_aspect("auto")
#ax.plot(time,v1)
plt.show()
