import numpy as np
import matplotlib.pyplot as plt

fname1="DAT_decay.out"
fp=open(fname1,"r")
xx=[]
amp=[]
for row in fp:
    dat=row.strip().split(",")
    xx.append(float(dat[0]))
    amp.append(float(dat[1]))

xx=np.array(xx)
amp1=np.array(amp)/amp[0]
fp.close()

fname2="Gss_decay.out"
fp=open(fname2,"r")
xx=[]
amp=[]
for row in fp:
    dat=row.strip().split(",")
    xx.append(float(dat[0]))
    amp.append(float(dat[1]))
xx=np.array(xx)
amp2=np.array(amp)/amp[0]
fp.close()

fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(xx,amp2,"r",linewidth=3)
ax.plot(xx,amp1,"b",linewidth=3)
ax.grid(True)
ax.tick_params(labelsize=14)
ax.set_xlim([0,43])
ax.set_ylim([0.2,1.2])
fig.savefig("decay.png",bbox_inches="tight")
plt.show()

