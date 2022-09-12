import numpy as np
import matplotlib.pyplot as plt
import sys
from PyQt5.QtWidgets import (QMainWindow, QFileDialog, QApplication)
#from PyQt5.QtWidgets import (QMainWindow, QTextEdit, QAction, QFileDialog, QApplication)
#from PyQt5.QtGui import QIcon

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
        fp.close()
    def show(self,ax,vmin=-0.2, vmax=0.2, cmap="jet",interpolation="none"):
        Xb=self.Xb;
        Xa=self.Xa;
        ext=[Xa[0],Xb[0],Xa[1],Xb[1]]
        im=ax.imshow(self.Z,aspect=1.0,extent=ext,cmap=cmap,origin="lower",vmin=vmin,vmax=vmax, interpolation=interpolation)
        ax.grid(True)
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
