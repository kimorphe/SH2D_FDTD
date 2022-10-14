#! /home/kazushi/anaconda3/bin/python
import sys
from PyQt5.QtWidgets import (
    QApplication, 
    QMainWindow,
    QWidget,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QHBoxLayout,
    QGridLayout,
    QFileDialog,
    QTextEdit,
    QLineEdit,
    QComboBox,
)
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from mkpng import FLD

class ImgPrms:
    def __init__(self):
        self.nfile=0
        self.fnames=[]
        self.cmap="jet"
        self.intpl="bilinear"
        self.vmin=0.0
        self.vmax=1.0
    def show(self):
        print("-------------------------")
        print("* number of files slected= ",self.nfile)
        print("* colormap= ",self.cmap)
        print("* interpolationp= ",self.intpl)
        print("* files: ")
        for fn in self.fnames:
            print(fn)
        print("-------------------------")

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        self.setWindowTitle("PNG Snapshot Generator")

        layout=QVBoxLayout()
        layout1=QHBoxLayout()
        layout2=QHBoxLayout()
        layout3=QHBoxLayout()

        btn1=QPushButton("Select Files")
        btn2=QPushButton("RUN")
        btn3=QPushButton("Apply")

        layout1.addWidget(btn1)
        layout1.addWidget(btn2)

        ledit1=QLineEdit("-0.15")
        ledit2=QLineEdit("0.15")
        lbl1=QLabel("min.")
        lbl2=QLabel("max.")
        layout2.addWidget(lbl1)
        layout2.addWidget(ledit1)
        layout2.addWidget(lbl2)
        layout2.addWidget(ledit2)


        cbox1=QComboBox()
        cbox2=QComboBox()
        cbox1.addItems(["jet","gray","gist_rainbow","ocean","nipy_spectral"])
        cbox2.addItems(["bilinear","bicubic","gaussian","spline16","none"])
        layout3.addWidget(cbox1)
        layout3.addWidget(cbox2)

        prms=ImgPrms()
        prms.vmin=float(ledit1.text())
        prms.vmax=float(ledit2.text())
        prms.cmap=cbox1.currentText()
        prms.itpl=cbox2.currentText()

        layout.addLayout(layout1)

        edit=QTextEdit()
        edit.setReadOnly(True)
        layout.addWidget(edit)

        layout2.addWidget(btn3)
        layout.addLayout(layout2)
        layout.addLayout(layout3)

        cntr_wdgt=QWidget()
        cntr_wdgt.setLayout(layout)

        self.setCentralWidget(cntr_wdgt)

        # Slot Assignment
        btn1.clicked.connect(self.btn1_pressed)
        btn3.clicked.connect(self.btn3_pressed)
        cbox1.currentTextChanged.connect(self.cbox1_changed)
        cbox2.currentTextChanged.connect(self.cbox2_changed)
        btn2.clicked.connect(self.RUN)
        #edit.clear()

        self.btn1=btn1
        self.btn2=btn2
        self.btn3=btn3
        self.edit=edit
        self.ledit1=ledit1
        self.ledit2=ledit2
        self.prms=prms

    def btn1_pressed(self,sig):
        fnames,typ=QFileDialog.getOpenFileNames(self,"Select files","./","*.out")
        for fn in fnames:
            self.edit.append(fn)
        nfile=len(fnames)
        self.edit.append(" -->Number of files selected ="+str(nfile))
        self.edit.append("---------------------------------")
        self.prms.fnames=fnames
        self.nfiles=len(fnames)

    def btn3_pressed(self,sig):  # Slot for "Apply" button 
        self.prms.vmin=float(self.ledit1.text())
        self.prms.vmax=float(self.ledit2.text())
        self.prms.show()

    def cbox1_changed(self,s):
        self.prms.cmap=s
    def cbox2_changed(self,s):
        self.prms.intpl=s

    def RUN(self):
        prms=self.prms
        v1=prms.vmin
        v2=prms.vmax
        cmap=prms.cmap
        intpl=prms.intpl

        if len(prms.fnames)==0:
            self.edit.append(" Files not selected")
            return()

        self.edit.append(" Generating PNG files.... ")
        fig=plt.figure()
        ax=fig.add_subplot(111)
        #ax.grid(True)
        v3=FLD()
        k=0
        for fn in prms.fnames:
            print(fn)
            v3.load(fn)
            if k==0:
                im=v3.show(ax,vmin=v1,vmax=v2, cmap=cmap, interpolation=intpl)
                divider=make_axes_locatable(ax)
                cax=divider.append_axes("right",size="5%",pad=0.05)
                plt.colorbar(im,cax=cax)
                ax.set_xlabel(r"$\it{x}$[mm]",fontsize=14)
                ax.set_ylabel(r"$\it{y}$[mm]",fontsize=14)
            else:
                im.set_data(v3.Z)
            stime=r"$t$={:6.2f}[micro sec]".format(v3.time)
            ax.set_title(stime,loc="center")
            fout=fn.replace("out","png")
            fig.savefig(fout,bbox_inches="tight")
            k+=1
        self.edit.append("... Done")


if __name__=="__main__":
    plt.rcParams["font.family"]="serif"
    plt.rcParams["font.serif"]=["Time New Roman"]+plt.rcParams["font.serif"]
    plt.rcParams["mathtext.fontset"]="stix"
    plt.rcParams["font.size"]=12
    #plt.rcParams["text.usetex"]=False
    app=QApplication(sys.argv)
    win=MainWindow()
    win.show()
    app.exec()
