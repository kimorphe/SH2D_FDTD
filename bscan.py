import numpy as np
import matplotlib.pyplot as plt


class WVs:
    def load(self,fname):
        fp=open(fname,"r")

        fp.readline()
        nele=int(fp.readline())

        fp.readline()
        dat=fp.readline().strip().split(",")
        Nt=int(dat[0])
        dt=float(dat[1])

        ndat=Nt*nele;
        amp=np.zeros(ndat)
        m=0;
        xrec=[]; yrec=[]
        for k in range(nele):
            dat=fp.readline()
            dat=dat.strip().split(",")
            xrec.append(float(dat[1]))
            yrec.append(float(dat[2]))
            for l in range(Nt):
                amp[m]=float(fp.readline())
                m+=1

        self.time=np.arange(Nt)*dt
        self.dt=dt
        self.amp=np.reshape(amp,[nele,Nt])
        self.Nt=Nt
        self.nele=nele
        self.xrec=xrec
        self.yrec=yrec
    def show_prms(self):
        print("nele=",self.nele)
        print("Nt,dt=",self.Nt,", ",self.dt)
        for k in range(self.nele):
            print(self.xrec[k],", ",self.yrec[k])

    def Aplot(self, ax):
        for k in range(self.nele):
            ax.plot(self.time, self.amp[k,:])

if __name__=="__main__":
    bwv=WVs()
    bwv.load("ary.out");

    fig=plt.figure()
    ax=fig.add_subplot(111)
    bwv.Aplot(ax)
    ax.grid(True)
    bwv.show_prms()

    plt.show()

