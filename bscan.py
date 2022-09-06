import numpy as np
import matplotlib.pyplot as plt


class WVs:
    def load(self,fname):
        fp=open(fname,"r")
        print(fp.readline())
        dat=fp.readline().strip().split(",")
        nele=int(dat[0]);
        Nt=int(dat[1])
        print(nele,Nt)

        ndat=Nt*nele;
        time=np.zeros(ndat)
        amp=np.zeros(ndat)
        m=0;
        for k in range(nele):
            print(fp.readline())
            for l in range(Nt):
                dat=fp.readline().strip().split(",")
                time[m]=float(dat[0])
                amp[m]=float(dat[1])
                m+=1
        self.time=np.reshape(time,[nele,Nt])
        self.amp=np.reshape(amp,[nele,Nt])
        self.Nt=Nt
        self.nele=nele
    def Aplot(self, ax):
        for k in range(self.nele):
            ax.plot(self.time[k,:], self.amp[k,:])

if __name__=="__main__":
    bwv=WVs()
    #bwv.load("ary0.out");
    bwv.load("ary1.out");

    fig=plt.figure()
    ax=fig.add_subplot(111)
    bwv.Aplot(ax)
    ax.grid(True)

    plt.show()

