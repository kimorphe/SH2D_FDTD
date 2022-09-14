import sys


from PyQt5 import QtCore, QtWidgets
from PyQt5.QtWidgets import (
        QApplication, 
        QMainWindow, 
        QPushButton,QWidget, 
        QVBoxLayout, 
        QHBoxLayout,
        QLabel,
        QSpinBox,
)
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg, NavigationToolbar2QT as NavigationToolbar
import numpy as np
import matplotlib


matplotlib.use("Qt5Agg")

class MplCanvas(FigureCanvasQTAgg):
    def __init__(self,parent=None, width=5, height=4, dpi=100):
        fig=plt.figure(figsize=(width,height),dpi=dpi)
        self.axes=fig.add_subplot(111)
        super(MplCanvas, self).__init__(fig)
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        #btn=QPushButton("button")
        sc=MplCanvas(self,width=5,height=4, dpi=100)
        x=np.linspace(0,2.*np.pi)
        y=np.sin(x)
        self.grph=sc.axes.plot(x,y)

        toolbar=NavigationToolbar(sc,self)

        layout=QVBoxLayout()
        layout2=QHBoxLayout()

        layout.addWidget(toolbar)
        layout.addWidget(sc)

        lbl1=QLabel("T")
        lbl2=QLabel("R")
        sboxT=QSpinBox()
        sboxT.setPrefix("T=")
        sboxT.setSingleStep(1)
        nsrc=4
        sboxT.setRange(1,nsrc)

        nrec=6
        sboxR=QSpinBox()
        sboxR.setSingleStep(1)
        sboxR.setPrefix("R=")
        sboxR.setRange(1,nrec)
        print(sboxR.value())

        sboxR.valueChanged.connect(self.val_changed)
        sboxT.valueChanged.connect(self.val_changed)
        #layout2.addWidget(lbl1)
        layout2.addWidget(sboxT)
        #layout2.addWidget(lbl2)
        layout2.addWidget(sboxR)
        layout.addLayout(layout2)

        self.sboxR=sboxR
        self.sboxT=sboxT

        wdgt=QWidget()
        wdgt.setLayout(layout)
        self.setCentralWidget(wdgt)
        self.x=x
        self.y=y
        self.sc=sc
        self.k=0

        Tp=10
        self.timer=QtCore.QTimer()
        self.timer.setInterval(Tp)
        #self.timer.timeout.connect(self.update)
        self.timer.timeout.connect(self.update2)
        self.timer.start()
    def val_changed(self):
        print("Spinbox value(T,R)=",self.sboxT.value(),self.sboxR.value())

    def update(self):
        dx=np.pi/20;
        k=self.k
        y=np.sin(self.x+dx*k)
        self.sc.axes.cla()  # clear the canvas
        self.sc.axes.plot(self.x,y)
        self.sc.draw()
        self.k+=1
    def update2(self):
        dx=np.pi/20;
        k=self.k
        y=np.sin(self.x+dx*k)
        self.grph[0].set_ydata(y)
        self.k+=1
        self.sc.draw()
        

if __name__=="__main__":
    app=QApplication(sys.argv)
    win=MainWindow()
    win.show()

    app.exec()

