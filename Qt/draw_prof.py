import matplotlib.pyplot as plt
import numpy as np

class BND:
    def __init__(self,fname,typ):
        fp=open(fname,"r")

        x=[];y=[];
        for row in fp:
            dat=row.strip().split(",")
            x.append(float(dat[0]))
            y.append(float(dat[1]))
        fp.close()

        self.x=np.array(x)
        self.y=np.array(y)
        self.typ=typ    # 0 for q1 boundary, 1 for q2 bounary

    def plot(self,ax):
        ax.plot(self.x,self.y,".")
    def plot2(self,ax,clr="k"):

        dx=0.15; 
        dy=0.15; 
        if self.typ==0:
            dx=0.
        if self.typ==1:
            dy=0.
        dx2=dx*0.5
        dy2=dy*0.5;

        nseg=len(self.x)
        for k in range(nseg):
            xc=self.x[k]
            yc=self.y[k]
            xx=[xc-dx2,xc+dx2]
            yy=[yc-dy2,yc+dy2]
            ax.plot(xx,yy,"-"+clr,linewidth=1.0)


if __name__=="__main__":
    bnd1=BND("q1bnd.dat",0)
    bnd2=BND("q2bnd.dat",1)


    fig=plt.figure()
    ax=fig.add_subplot(111)
    bnd1.plot(ax,clr="m")
    bnd2.plot(ax,clr="g")
    #bnd1.plot(ax)
    #bnd2.plot(ax)
    ax.grid(True)
    ax.set_aspect(1.0)
    plt.show()



