import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

import sys
from PyQt5.QtWidgets import (QMainWindow, QFileDialog, QApplication)



class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        #sys.exit()
    def showDialog(self):
        fnames,typ=QFileDialog.getOpenFileNames(self,"Select files","./","*.out")
        #fname,=QFileDialog.getOpenFileName(self,"Open file","/home")
        return(fnames)

class FLD:
    def load(self,fname):
        fp=open(fname,"r")

        cmt=fp.readline()
        time=float(fp.readline())
        cmt=fp.readline()
        dat=fp.readline().strip().split(",");
        x=float(dat[0]); y=float(dat[1])
        Xa=[x,y]

        cmt=fp.readline()
        dat=fp.readline().strip().split(",");
        x=float(dat[0]); y=float(dat[1])
        Xb=[x,y]

        cmt=fp.readline()
        dat=fp.readline().strip().split(",");
        Nx=int(dat[0])
        Ny=int(dat[1])

        cmt=fp.readline()
        Z=fp.readlines()
        Z=np.array(Z)
        Z=Z.astype(float)
        Z=np.reshape(Z,[Nx,Ny])

        self.Nx=Nx; self.Ny=Ny
        self.Z=np.transpose(Z)
        self.Xa=np.array(Xa)
        self.Xb=np.array(Xb)
        self.time=time
        fp.close()
    def show(self,ax,vmin=-0.2, vmax=0.2, cmap="jet",interpolation="none"):
        Xb=self.Xb;
        Xa=self.Xa;
        ext=[Xa[0],Xb[0],Xa[1],Xb[1]]
        im=ax.imshow(self.Z,aspect=1.0,extent=ext,cmap=cmap,origin="lower",vmin=vmin,vmax=vmax, interpolation=interpolation)
        #ax.grid(True)
        ax.set_xlim([Xa[0],Xb[0]])
        ax.set_ylim([Xa[1]-10,Xb[1]+10])
        ax.set_facecolor("k")
        return(im)

if __name__=="__main__":

    app=QApplication(sys.argv)
    window=MainWindow()
    #window.show()
    #app.exec()

    fnames=window.showDialog()
    print("----- Selected Files ------")
    k=1
    for fn in fnames:
        print(str(k)+" "+fn)
        k+=1
    print("---------------------------")
    #sys.exit(app.exec_())

    fig=plt.figure()
    ax=fig.add_subplot(111)

    F=FLD()
    k=0
    for fn in fnames:
        F.load(fn)
        if k==0:
            im=F.show(ax)
        else:
            im.set_data(F.self.Z)
        fn=fn.replace("out","png")
        fig.savefig(fn,bbox_inches="tight")
        print(fn)
