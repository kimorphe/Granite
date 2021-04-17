import numpy as np
import matplotlib.pyplot as plt


freq=np.linspace(0,2,301)

fig=plt.figure()
ax1=fig.add_subplot(311)
ax2=fig.add_subplot(312)
ax3=fig.add_subplot(313)

rho=1.0
E1=1
E2=2


Z1=E1
Z2=E2
omg=np.pi*2*freq;
eta=1
Z3=(-1j*omg)*eta


# 
M2=Z2+Z3 # Voigt
M1=Z1*Z3/(Z1+Z3) # Maxwell

etas=[0.05,0.1,0.5,1]
for eta in etas:
    Z3=(-1j*omg)*eta
    M1=Z1*Z3/(Z1+Z3) # Maxwell
    M2=Z2+Z3 # Voigt
    M3=Z1*M2/(Z1+M2) # Kelvin

    kw1=np.sqrt(rho/M1)*omg
    kw2=np.sqrt(rho/M2)*omg
    kw3=np.sqrt(rho/M3)*omg
    ax1.plot(freq,np.real(kw1),label=str(eta))
    ax2.plot(freq,np.real(kw2),label=str(eta))
    ax3.plot(freq,np.real(kw3),label=str(eta))
ax1.grid(True); ax1.legend()
ax2.grid(True); ax2.legend()
ax3.grid(True); ax3.legend()
ax1.set_title("Maxwell")
ax2.set_title("Voigt")
ax3.set_title("Kelvin")


plt.show()
