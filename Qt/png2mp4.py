import sys
import os
from PyQt5 import QtCore
from PyQt5.QtWidgets import(
        QApplication,
        QMainWindow,
        QWidget,
        QHBoxLayout,
        QVBoxLayout,
        QLabel,
        QPushButton,
        QSpinBox,
        QFileDialog,
        QTextEdit,
        QLineEdit,
)

import subprocess

class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("png -->  mpeg converter")

        Vlayout=QVBoxLayout()
        Hlayout=QHBoxLayout()
        btn1=QPushButton("Select Directory")
        btn1.clicked.connect(self.select_directory)

        dsply=QTextEdit()
        dsply.setReadOnly(True)
        Vlayout.addWidget(dsply)

        lbl=QLabel("Frame Rate")
        ledt=QLineEdit("10")
        ledt.setAlignment(QtCore.Qt.AlignCenter)

        lbl2=QLabel("Write to=")
        ledt2=QLineEdit("v3")
        ledt2.setAlignment(QtCore.Qt.AlignRight)
        lbl3=QLabel(".mp4")

        btn2=QPushButton("Generate mp4 Movie ")
        btn2.setEnabled(False)
        btn2.clicked.connect(self.run_ffmpeg)

        Hlayout.addWidget(btn1)
        Hlayout.addWidget(lbl)
        Hlayout.addWidget(ledt)

        Hlayout.addWidget(lbl2)
        Hlayout.addWidget(ledt2)
        Hlayout.addWidget(lbl3)

        #Hlayout.addWidget(btn2)

        Vlayout.addLayout(Hlayout)
        Vlayout.addWidget(btn2)
        cntrl_widget=QWidget()
        cntrl_widget.setLayout(Vlayout)
        self.setCentralWidget(cntrl_widget)


        self.btn2=btn2  #execute button
        self.dsply=dsply
        self.ledt=ledt
        self.ledt2=ledt2

    def run_ffmpeg(self):
        path=self.dir_name+"/v3_%d.png"
        frt=self.ledt.text()
        fn=self.ledt2.text()

        fout=self.dir_name+"/"+fn+".mp4"
        #cmd=["ffmpeg","-r",frt,"-i",path,"-vcodec","h264",fout]
        opt1=["-vcodec", "h264", "-y"]
        opt2=["-pix_fmt", "yuv420p", "-vf", "scale=trunc(iw/2)*2:trunc(ih/2)*2"]
        cmd=["ffmpeg","-r",frt,"-i",path]+opt1+opt2+[fout]
        print(cmd)
        #ffmpeg -r 5 -i v3_%d.png -vcodec h264 -pix_fmt yuv420p -vf "scale=trunc(iw/2)*2:trunc(ih/2)*2" out.mp4

        self.dsply.append("Running ffmpeg ...\n")
        self.dsply.append(" Frame rate= "+frt)
        subprocess.run(cmd,text=True)
        self.dsply.append(" MP4 movie "+fout+" created")

    def select_directory(self):
        dir_name=QFileDialog.getExistingDirectory(self,"Select Directory","./",QFileDialog.ShowDirsOnly)
        self.dir_name=dir_name
        self.btn2.setEnabled(True)

        self.dsply.clear()
        self.dsply.append("## Current data folder:")
        self.dsply.append("  -->"+dir_name)
        """
        nfile=0
        head="v3_";
        tail=".out"
        while os.path.exists(head+str(nfile)+tail):
            nfile+=1
        self.nfile=nfile
        """

if __name__=="__main__":
    app=QApplication(sys.argv)
    win=MainWindow()
    win.show()
    app.exec()

