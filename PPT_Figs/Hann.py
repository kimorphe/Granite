import numpy as np
import matplotlib.pyplot as plt

x=np.linspace(0,6,201)

s=1.0

xi=x/s-1
W=0.5*(1.-np.cos(np.pi*(xi+1)))

indx=np.argwhere(xi>1)
W[indx]=0.0

fig=plt.figure()
ax=fig.add_subplot(111)

ax.plot(xi,W)
ax.grid(True)

plt.show()
