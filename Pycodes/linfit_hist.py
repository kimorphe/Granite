import numpy as np
import matplotlib.pyplot as plt

fp=open("linfit_K.out","r");
a=[]    # gradient
b=[]    # intercept  
c=[]    # mean cp
k=0;
for row in fp:
    dat=row.strip().split(",")
    a.append(float(dat[0]));
    b.append(float(dat[1]));
    c.append(float(dat[2]));

a=np.array(a)
b=np.array(b)
c=np.array(c)

fp.close()
fig=plt.figure()
ax=fig.add_subplot(121)
bx=fig.add_subplot(122)

ax.hist(a,range=[-3,3],bins=30)
bx.hist(c,range=[3,8],bins=30)
ax.grid(True)
bx.grid(True)

plt.show()
