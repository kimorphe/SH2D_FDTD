import matplotlib.pyplot as plt
import numpy as np

class ARRAY:
    def __init__(self,nele):
        self.nele=nele
        self.tdly=np.zeros(nele)
        self.amp=np.ones(nele)
        self.nmeas=0
        self.fout=""
    def deploy(self,pitch,width):
        W=(self.nele-1)*pitch
        self.xcod=np.arange(self.nele)*pitch-0.5*W
        self.ycod=np.zeros(self.nele)
        self.width=width
    def shift(self,x0,y0):
        self.xcod+=x0
        self.ycod+=y0
    def focus(self,xf,phase_vel):
        rx=xf[0]-self.xcod
        ry=xf[1]-self.ycod
        tof=np.sqrt(rx*rx+ry*ry)/phase_vel
        tmax=np.max(tof)
        self.tdly=tmax-tof

        self.fout+="# active, A0, delay\n"
        actv=1
        for k in range(self.nele):
            self.fout+="{},{},{}\n".format(actv,self.amp[k],self.tdly[k])
        self.nmeas+=1


    def steer(self,th_in_deg,phase_vel):
        th=np.deg2rad(th_in_deg)
        if th>0.0:
            self.tdly=np.sin(th)*(self.xcod-self.xcod[0])/phase_vel
        else:
            self.tdly=np.sin(th)*(self.xcod-self.xcod[-1])/phase_vel

        self.fout+="# active, A0, delay\n"
        actv=1
        for k in range(self.nele):
            self.fout+="{},{},{}\n".format(actv,self.amp[k],self.tdly[k])
        self.nmeas+=1

    def fwrite_elem(self):
        fout=""
        cmt="# nsrc(number of source element)"
        fout+=cmt+"\n"
        fout+=str(nele)+"\n"
        nml=1
        wvID=0
        th_in=0.0   # intraelement delay
        cmt="# type, xy, wd, nml, wvID, th_in"
        cmt+="(type=1:v1, 2:v2, posi., width, normal, wvfm No, th_in[deg])"
        fout+=cmt+"\n"
        for k in range(self.nele):
            fout+="{},{},{},{},{},{}\n".format(2,self.xcod[k],self.width,nml,wvID,th_in)
        print(fout)

        fp=open("ary_elem.txt","w")
        fp.write(fout)
        fp.close()

    def fwrite_delay(self):
        fout=""
        cmt="# nsrc(number of source element)"
        fout+=cmt+"\n"
        fout+=str(nele)+"\n"
        fout+="# N_meas (Array)\n"
        fout+=str(self.nmeas)+"\n"

        #fout+="# active, A0, delay\n"
        #actv=1
        #for k in range(self.nele):
        #    fout+="{},{},{}\n".format(actv,self.amp[k],self.tdly[k])

        fout+=self.fout
        print(fout)
        fp=open("ary_delay.txt","w")
        fp.write(fout)
        fp.close()


    

if __name__=="__main__":

    xf=[10., 20.]
    ycod=0.0

    nele=32
    ary=ARRAY(nele)
    pitch=1.0
    width=0.6
    ary.deploy(pitch,width)

    vel=3.0
    ary.focus(xf,vel)
    print(ary.tdly)
    ary.steer(-60.,vel)
    print(ary.tdly)
    ary.fwrite_elem()
    ary.fwrite_delay()

    """
    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.plot(xcod,tof,"o-")
    ax.plot(xcod,tdly,"s-")
    ax.plot(xcod,tof+tdly,"-k")
    ax.grid(True)
    plt.show()
    """
