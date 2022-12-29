#! /home/kazushi/anaconda3/bin/python
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

fname="saft.out";


fp=open(fname,"r");

fp.readline(); # computational domain
Xa=list(map(float,fp.readline().lstrip().split(" ")));

fp.readline(); # computational domain
Xb=list(map(float,fp.readline().lstrip().split(" ")));

fp.readline(); # computational domain
Ndiv=list(map(int,fp.readline().lstrip().split(" ")));

fp.readline();	# Imaging area
K=list(map(float,fp.readlines()));	# Imaging area
K=np.transpose(np.reshape(K,Ndiv))

print(np.shape(K))


fig=plt.figure()
ax=fig.add_subplot(111)
im=ax.imshow(K,origin="lower",extent=[Xa[0],Xb[0],Xa[1],Xb[1]],interpolation="bilinear",cmap="jet")

plt.show()
