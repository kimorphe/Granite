import numpy as np

import matplotlib.pyplot as plt

fname="inwv.out"
fp=open(fname,"r")

tt=[]
amp=[]
for row in fp:
    dat=row.strip().split()
    tt.append(float(dat[0]))
    amp.append(float(dat[1]))

tt=np.array(tt)
amp=np.array(amp)

fig=plt.figure()
ax=fig.add_subplot(111)
ax.plot(tt,amp,"b",linewidth=3)
ax.grid(True)
ax.tick_params(labelsize=16)
ax.set_xlim([0,4])
fig.savefig("inwv.png",bbox_inches="tight")
plt.show()

