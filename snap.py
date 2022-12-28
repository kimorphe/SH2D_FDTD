import numpy as np
import matplotlib.pyplot as plt


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
        #wx=float(dat[0]); wy=float(dat[1])
        #Wd=[wx,wy]
        x=float(dat[0]); y=float(dat[1])
        Xb=[x,y]

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
        self.Z=np.transpose(Z)
        self.Xa=np.array(Xa)
        self.Xb=np.array(Xb)
        #self.Wd=np.array(Wd)
    def show(self,ax):
        Xb=self.Xb;
        Xa=self.Xa;
        ext=[Xa[0],Xb[0],Xa[1],Xb[1]]
        ax.imshow(self.Z,aspect=1.0,extent=ext,cmap="jet",origin="lower",vmin=-0.2,vmax=0.2)
        ax.grid(True)

if __name__=="__main__":
    F=FLD()


    fig=plt.figure()
    ax=fig.add_subplot(111)

    nums=np.arange(0,61,2)
    print(nums)
    for k in nums:
        fname="v3_"+str(k)+".out"
        F.load(fname)
        F.show(ax)
        fname=fname.replace("out","png")
        fig.savefig(fname,bbox_inches="tight")

    plt.show()

