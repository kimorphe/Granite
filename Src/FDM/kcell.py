#! /home/kazushi/Enthought/Canopy_64bit/User/bin/python
import numpy as np
import matplotlib.pyplot as plt


fname="kcell.dat"; 

fp=open(fname,"r");

fp.readline(); # computational domain
xa=list(map(float,fp.readline().lstrip().split(" ")));
xb=list(map(float,fp.readline().lstrip().split(" ")));

fp.readline();	# Physical domain
Xa=list(map(float,fp.readline().lstrip().split(" ")));
Xb=list(map(float,fp.readline().lstrip().split(" ")));

fp.readline();	# Imaging area
Ya=list(map(float,fp.readline().lstrip().split(" ")));
Yb=list(map(float,fp.readline().lstrip().split(" ")));

fp.readline();	# Imaging area
Ndiv=list(map(int,fp.readline().lstrip().split(" ")));
Ng=Ndiv[0]*Ndiv[1];

fp.readline();	# Imaging area
K=list(map(float,fp.readlines()));	# Imaging area
K=np.transpose(np.reshape(K,Ndiv))

print(np.shape(K));

fig=plt.figure();
ax=fig.add_subplot(111)
ax.imshow(K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],cmap="jet",interpolation="none");
plt.show();






fp.close();
