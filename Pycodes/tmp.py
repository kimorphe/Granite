import numpy as np
import matplotlib.pyplot as plt

fig=plt.figure()
ax=fig.add_subplot(111)

x=np.linspace(0,100,2000)
tb=50
wt=2
a=1
x90=10
a=np.log(0.9/(1-0.9))/x90
y=1/(1+np.exp(-a*(x-tb)))
ax.plot(x,y)
ax.grid(True)

plt.show()
