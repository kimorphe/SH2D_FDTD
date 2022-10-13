import numpy as np
import matplotlib.pyplot as plt
import kcell


class TRdat:
    def load(self,fname):
        fp=open(fname)

        NL=6;
        for k in range(NL):
            row=fp.readline()
            if k==4:
                ng=row.replace("# ng=","")
                ng=int(ng)
                print("ng=",ng)

        xe=[]
        ye=[]
        for k in range(ng):
            dat=fp.readline()
            dat=dat.strip().split(",")
            xe.append(float(dat[2]))
            ye.append(float(dat[3]))
        fp.close()

        self.xe=np.array(xe)
        self.ye=np.array(ye)
    def draw(self,ax):
        ax.plot(self.xe,self.ye,"y",linewidth=3)

if __name__=="__main__":

    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times New Roman'] + plt.rcParams['font.serif']
    plt.rcParams['text.usetex']==False
    plt.rcParams['mathtext.fontset']='stix'
    plt.rcParams['font.size']=12

    tr=TRdat()
    tr.load("tr_elems.out")

    K=kcell.KCELL()
    fname="kcell.dat"
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.set_facecolor("k")
    fsz=14
    #ax.tick_params(labelsize=fsz)
    ax.set_xlabel(r'$x$[mm]',fontsize=fsz)
    ax.set_ylabel("$y$[mm]",fontsize=fsz)
    K.load(fname)
    K.show(ax)

    tr.draw(ax)
    ax.set_ylim([-5,50])
    

    fname="model.png"
    fig.savefig(fname,bbox_inches="tight")

    plt.show()


