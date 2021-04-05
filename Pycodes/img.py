from PIL import Image
import numpy as np
import matplotlib.pyplot as plt

im=np.array(Image.open("core_top_face.png"))
print(type(im))
print(im.dtype)
print(im.shape)

fig=plt.figure()
ax=fig.add_subplot(111)
ax.imshow(im)

Xa=[300,300]
Xb=[900,900]

xlim=ax.get_xlim()
ylim=ax.get_xlim()

ax.vlines(Xa[0],ylim[0],ylim[1])
ax.vlines(Xb[0],ylim[0],ylim[1])
ax.hlines(Xa[1],xlim[0],xlim[1])
ax.hlines(Xb[1],xlim[0],xlim[1])

IM=im[Xa[0]:Xb[0],Xa[1]:Xb[1],:]
print(np.shape(IM))

fig2=plt.figure()
bx=fig2.add_subplot(121)
bx.imshow(IM[:,:,:])
cx=fig2.add_subplot(122)
ix=325
Rx=IM[ix,:,0]
Gx=IM[ix,:,1]
Bx=IM[ix,:,2]

cx.plot(IM[ix,:,0],"-r",markersize=1,alpha=0.5)
cx.plot(IM[ix,:,1],"-g",markersize=1,alpha=0.5)
cx.plot(IM[ix,:,2],"-b",markersize=1,alpha=0.5)
#cx.plot(IM[:,:,0],np.mean(IM,axis=2),"o",markersize=1,alpha=0.3)
#cx.plot(IM[:,:,0],IM[:,:,1],".b",markersize=1,alpha=0.3)
xlim=bx.get_xlim()
ylim=bx.get_ylim()
bx.hlines(ix,xlim[0],xlim[1])
cx.grid(True)


fig3=plt.figure()
dx1=fig3.add_subplot(411)
dx2=fig3.add_subplot(412)
dx3=fig3.add_subplot(413)
dx4=fig3.add_subplot(414)
ndat=np.shape(IM)
R=np.reshape(IM[:,:,0],[ndat[0]*ndat[1]])
G=np.reshape(IM[:,:,1],[ndat[0]*ndat[1]])
B=np.reshape(IM[:,:,2],[ndat[0]*ndat[1]])
dx1.hist(R,bins=40,color="r")
dx2.hist(G,bins=40,color="b",alpha=0.5)
dx3.hist(B,bins=40,color="g",alpha=0.75)
dx1.grid(True)
dx2.grid(True)
dx3.grid(True)

K=[170,155,155]
Na=[150,175,200]
Q=[100,110,120]
Bt=[50,50,50]

K=[150,155,167]
Na=[131,137,148]
Q=[89,97,106]
Bt=[62,68,76]

D=np.reshape([K,Na,Q,Bt],[1,4,3])

dx1.set_xlim([0,256])
dx2.set_xlim([0,256])
dx3.set_xlim([0,256])
dx4.imshow(D)


fp=open("rgb.dat","w")
print(len(R))
fp.write(str(len(R))+"\n");
for k in range(len(R)):
    dat=str(R[k])+", "+str(G[k])+", "+str(B[k])+"\n"
    fp.write(dat)
fp.close()

fp=open("rgb1d.dat","w")
n=len(Rx)
fp.write(str(n)+"\n");
for k in range(n):
    dat=str(Rx[k])+", "+str(Gx[k])+", "+str(Bx[k])+"\n"
    fp.write(dat)
fp.close()

plt.show()
