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

        self.time=np.arange(Nt)*dt  #time axis
        self.dt=dt  # time increment
        self.amp=np.reshape(amp,[nele,Nt])  # wave data
        self.Nt=Nt  # number of time steps
        self.nele=nele  # number of array elements
        self.xrec=xrec  # element position (x)
        self.yrec=yrec  # element position (y)
    def show_prms(self):
        print("nele=",self.nele)
        print("Nt,dt=",self.Nt,", ",self.dt)
        for k in range(self.nele):
            print(self.xrec[k],", ",self.yrec[k])
    def Aplot(self, ax,nums):
        for k in nums:
            ax.plot(self.time, self.amp[k,:],label="R="+str(k))
        ax.set_xlim(self.time[0],self.time[-1])
    def Bplot(self, ax,v1=-0.05,v2=0.05):
        ext=[self.time[0],self.time[-1],self.xrec[0],self.xrec[-1]]
        ax.imshow(self.amp,extent=ext,cmap="jet",aspect="auto",origin="lower",interpolation="none",vmin=v1,vmax=v2)

    def correlate(self,num,T0):
        aref=np.zeros(self.Nt)
        icut=int(T0*8/self.dt)
        aref[0:icut]+=self.amp[num,0:icut]
        #aref/=np.sqrt(np.sum(aref*aref))
        aref/=(np.sum(aref*aref))
        self.aref=aref

        indx=self.Nt-np.array(range(self.Nt))-1
        cor_wv=np.zeros(np.shape(self.amp))
        for k in range(self.nele):
            y=np.correlate(self.aref,self.amp[k,:],mode="full");
            cor_wv[k,:]=y[indx]
            #ax2.plot(self.time,cor_wv[k,:])
        self.cor_wv=cor_wv
    def plot_CorWv_B(self,ax,v1=-0.1,v2=0.1):
        ext=[self.time[0],self.time[-1],self.xrec[0],self.xrec[-1]]
        ax.imshow(self.cor_wv,extent=ext,cmap="jet",aspect="auto",origin="lower",interpolation="none",vmin=v1,vmax=v2)
    def plot_CorWv_A(self,ax,nums):
        for k in nums:
            ax.plot(self.time, self.cor_wv[k,:],label="R="+str(k))
        ax.set_xlim(self.time[0],self.time[-1])
    def export_CorWv(self,fname):
        TXT="#nele\n"
        TXT+=str(self.nele)+"\n"
        TXT+="# Nt, dt"+"\n"
        TXT+=str(self.Nt)+", "+str(self.dt)+"\n"
        for k in range(self.nele):
            x=str(self.xrec[k])
            y=str(self.yrec[k])
            TXT+=("# e="+str(k)+", "+x+", "+y+"\n")
            for j in range(self.Nt):
                TXT+=(str(self.cor_wv[k,j])+"\n")

        fp=open(fname,"w")
        fp.write(TXT)
        fp.close()


if __name__=="__main__":
    f0=1.0  # signal frequency
    T0=1./f0 # period

    bwv=WVs()   # B-scan class instance

    #T_scan=[5,15]
    T_scan=range(N_scan)
    N_scan=32
    for Tnum in T_scan:
        fname="T"+str(Tnum)+"/ary.out"
        fnout="T"+str(Tnum)+"/ary_cor.out"
        print(fname+" --> "+fnout)
        bwv.load(fname) # load B-scan data
        bwv.correlate(Tnum,T0) # define reference signal and get correlation functions
        #bwv.show_prms() # show B-scan properties

        # plot waveforms
        if(len(T_scan)<3):
                enums=[0,15,31,Tnum]    # A-scan plot 
                fig1=plt.figure()
                ax1=fig1.add_subplot(211)
                bx1=fig1.add_subplot(212)
                bwv.Aplot(ax1,enums)
                bwv.Bplot(bx1)
                ax1.grid(True)
                ax1.legend()
                ax1.set_title("T="+str(Tnum)+" (raw data)")
                fig2=plt.figure()
                ax2=fig2.add_subplot(211)
                bx2=fig2.add_subplot(212)
                bwv.plot_CorWv_A(ax2,enums)
                bwv.plot_CorWv_B(bx2)
                ax2.grid(True)
                ax2.legend()
                ax2.set_title("T="+str(Tnum)+" (correlated waveforms)")
                plt.show()

        bwv.export_CorWv(fnout)

