import sys
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import matplotlib
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
matplotlib.use("Qt5Agg")

from PyQt5.QtWidgets import (
        QApplication,
        QMainWindow,
        QWidget,
        QVBoxLayout,
        QHBoxLayout,
        QPushButton,
        QLabel,
        QFileDialog,
        QSpinBox,
        QRadioButton,
) 
class MplCanvas(FigureCanvasQTAgg):
    def __init__(self,parent=None, width=5, height=4, dpi=100):
        fig=plt.figure(figsize=(width,height),dpi=dpi)
        self.axes=fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)

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

    def Aplots(self, ax):
        for k in range(self.nele):
            ax.plot(self.time, self.amp[k,:])
    def Aplot(self, ax,num):
        ax.plot(self.time, self.amp[num,:])

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        btn1=QPushButton("select directory") 
        btn1.clicked.connect(self.select_dir)
        canvas=MplCanvas(self,width=5,height=4, dpi=100)
        toolbar=NavigationToolbar(canvas,self)

        layout=QVBoxLayout()
        layout.addWidget(btn1)
        layout.addWidget(toolbar)
        layout.addWidget(canvas)

        layout2=QHBoxLayout()
        spbT=QSpinBox(); 
        spbT.setRange(1,1); 
        spbT.setSingleStep(1); 
        spbT.setPrefix("T="); 
        spbT.valueChanged.connect(self.draw_Ascan)

        spbR=QSpinBox()
        spbR.setRange(1,1)
        spbR.setSingleStep(1)
        spbR.setPrefix("R=");
        spbR.valueChanged.connect(self.draw_Ascan)

        btn2=QPushButton("clear")
        btn2.clicked.connect(self.clear_canvas)
        layout2.addWidget(spbT)
        layout2.addWidget(spbR)

        rbtn=QRadioButton()
        lbl=QLabel("freez y-limit")
        #rbtn.setChecked(True)
        layout2.addWidget(lbl)
        layout2.addWidget(rbtn)
        layout2.addWidget(btn2)

        layout.addLayout(layout2)

        wgt=QWidget()
        wgt.setLayout(layout)
        self.setCentralWidget(wgt)

        self.canvas=canvas
        self.spbT=spbT
        self.spbR=spbR
        self.rbtn=rbtn

    def clear_canvas(self):
        self.canvas.axes.cla()
        self.canvas.draw()
    def draw_Ascan(self):
        self.canvas.axes.plot()
        T=self.spbT.value()-1 # get integer type value
        R=self.spbR.value()-1 # get integer type value
        print("draw Ascan called")
        bw=self.bwvs[T]
        ylim=self.canvas.axes.get_ylim()
        bw.Aplot(self.canvas.axes,R)
        #ylim=[-0.1,0.1]
        if self.rbtn.isChecked():
            self.canvas.axes.set_ylim(ylim)

        self.canvas.draw() 

    def select_dir(self):
        dir_name=QFileDialog.getExistingDirectory(self,"Select Directory","./",QFileDialog.ShowDirsOnly)
        print("Selected directory=",dir_name)
        self.dir_name=dir_name
        #print("Files & Directies: ", os.listdir(dir_name))

        fn_ary=dir_name+"/"+"array.inp"
        fn_src=dir_name+"/"+"src.inp"
        ary=ary_prms(fn_ary)
        src=src_prms(fn_src)
        ary.show()
        src.show()

        bwvs=[]
        for k in range(ary.N_meas):
            fname=dir_name+"/T"+str(k)+"/ary.out"
            print(fname)
            bwv=WVs()
            bwv.load(fname)
            bwv.show_prms()
            bwvs.append(bwv)
        #fig=plt.figure()
        #ax=fig.add_subplot(111)
        T=0
        R=1
        bw=bwvs[T]
        #bw.Aplot(self.canvas.axes)
        #self.canvas.draw()

        self.spbT.setRange(1,ary.N_meas)
        self.spbR.setRange(1,ary.nele)

        T=self.spbT.value()-1 # get integer type value
        R=self.spbR.value()-1 # get integer type value
        bw=bwvs[T]
        bw.Aplot(self.canvas.axes,R)
        self.canvas.draw()

        self.bwvs=bwvs


        """
        dirs=glob.glob(dir_name+"/*/")
        print("Directories:", dirs)
        for dat in dirs:
            print("Base name=",os.path.basename(dat.rstrip("/")))
            fname=os.path.basename(dat.rstrip("/"))+"/ary.out"
            fname=dir_name+"/"+fname
            print("File name",fname)
            
            bwv.load(fname)
            bwv.show_prms()
        """

class ary_prms:
    def __init__(self,fname):
        fp=open(fname,"r")
        fp.readline()
        nele=int(fp.readline())

        fp.readline()
        N_meas=int(fp.readline())

        actv=[]
        A0=[]
        tdly=[]
        for k in range(N_meas):
            fp.readline()
            for j in range(nele):
                dat=fp.readline()
                dat=dat.strip().split(",")
                actv.append(int(dat[0]))
                A0.append(float(dat[1]))
                tdly.append(float(dat[2]))
        fp.close()

        self.nele=nele
        self.N_meas=N_meas
        self.actv=actv
        self.A0=A0
        self.tdly=tdly
    def show(self):
        print("-------- Array Setting (array.inp)----------")
        print("nele=",self.nele)
        print("N_meas=",self.N_meas)
        l=0
        for k in range(self.N_meas):
            for j in range(self.nele):
                print("actv, A0, tdly=",self.actv[l],self.A0[l],self.tdly[l])
                l+=1

class src_prms:
    def __init__(self,fname):
        fp=open(fname,"r")
        fp.readline()
        nsrc=int(fp.readline())
        fp.readline()
        typ=[]
        xy=[]
        wd=[]
        nml=[]
        wvID=[]
        th_in=[]
        for k in range(nsrc):
            dat=fp.readline()
            dat=dat.strip().split(",")
            print(dat)
            typ.append(int(dat[0]))
            xy.append(float(dat[1]))
            wd.append(float(dat[2]))
            nml.append(int(dat[3]))
            wvID.append(int(dat[4]))
            th_in.append(float(dat[5]))
        fp.close()

        self.nsrc=nsrc
        self.typ=typ
        self.xy=xy
        self.wd=wd
        self.nml=nml
        self.wvID=wvID
        self.th_in=th_in

    def show(self):
        print("-------- Source Setting (src.inp) ----------")
        print("nsrc=",self.nsrc)
        print("typ=",self.typ)
        print("xy=",self.xy)
        print("wd=",self.wd)
        print("nml=",self.nml)
        print("wvID=",self.wvID)
        print("th_in=",self.th_in)

if __name__=="__main__":
    app=QApplication(sys.argv)
    win=MainWindow()
    win.show()
    app.exec()

