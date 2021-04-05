import numpy as np
import matplotlib.pyplot as plt


freq=np.linspace(0.5,1.5,1001);
omg=freq*2*np.pi
alpha=1.0
beta=1
c0=3.0

fig=plt.figure()
ax=fig.add_subplot(111)

Bs=[0.1,1,10]
As=[0.6,0.7,0.8,0.9,0.95,0.98,1.0]
#for beta in Bs:
for alpha in As:
    Z=np.sqrt(1+beta*((-1j*omg)**alpha))
    c=1/np.real(1/Z)

    M=1+beta*((-1j*omg)**alpha)
    Mb=np.abs(M)
    dlt=np.angle(M)
    c=np.sqrt(Mb)/np.cos(dlt*0.5)

    #ax.plot(freq,c*c0,label=str(alpha))
    #ax.loglog(omg,c,label=str(alpha))
    ax.plot(np.log10(omg),np.log10(c-0.5*c0),label=str(alpha));
    ax.set_aspect(1.0)

ax.grid(True)
ax.legend()

plt.show()
