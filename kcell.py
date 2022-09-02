import numpy as np
import matplotlib.pyplot as plt


class KCELL:
    def load(self,fname):
        fp=open(fname,"r")

        cmt=fp.readline()
        dat=fp.readline().strip().split(",");
        x=float(dat[0]); y=float(dat[1])
        Xa=np.array([x,y])

        dat=fp.readline().strip().split(",");
        x=float(dat[0]); y=float(dat[1])
        Xb=np.array([x,y])
        Wd=Xb-Xa;

        cmt=fp.readline()
        dat=fp.readline().strip().split(",");
        Nx=int(dat[0])
        Ny=int(dat[1])
        Z=[]
        cmt=fp.readline()
        for row in fp:
            Z.append(float(row))

        Z=np.array(Z)
        Z=np.reshape(Z,[Nx,Ny])

        self.Nx=Nx
        self.Ny=Ny
        self.Ndiv=np.array([Nx,Ny])
        self.Z=np.transpose(Z)
        self.Xa=Xa
        self.Xb=Xb
        self.Wd=Wd
        self.dx=self.Wd/self.Ndiv
    def show(self,ax):
        dx=self.dx;
        Xa=self.Xa;
        Xb=self.Xb;
        ext=[Xa[0],Xb[0],Xa[1],Xb[1]]
        ax.imshow(self.Z,aspect=1.0,extent=ext,cmap="jet",origin="lower")
        ax.grid(True)

if __name__=="__main__":
    K=KCELL()


    fname="kcell.dat"
    fig=plt.figure()
    ax=fig.add_subplot(111)
    K.load(fname)
    K.show(ax)
    #fname=fname.replace("out","png")
    #fig.savefig(fname,bbox_inches="tight")

    plt.show()

