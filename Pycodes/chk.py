import numpy as np
import matplotlib.pyplot as plt
import kw


fig=plt.figure()
ax1=fig.add_subplot(121)
ax2=fig.add_subplot(122)

awv=kw.AWV()
f0=1.0; Nt=501;Td=50
nbrst=4
awv.burst(nbrst,f0,Nt,Td)

awv.plot_ascan(ax1)
awv.plot_FFT(ax2)

T0=1/f0
tb=T0*nbrst*0.5
nsig=6
sig=4*T0/nsig
awv.Amod_Gauss(tb,sig)
awv.plot_ascan(ax1)
awv.plot_FFT(ax2)

plt.show()


