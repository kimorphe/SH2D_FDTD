import sys
import os
import glob

from PyQt5.QtWidgets import (
        QApplication,
        QMainWindow,
        QWidget,
        QVBoxLayout,
        QHBoxLayout,
        QPushButton,
        QLabel,
        QFileDialog,
) 


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()

        btn1=QPushButton("select directory") 
        btn1.clicked.connect(self.select_dir)

        layout=QVBoxLayout()
        layout.addWidget(btn1)


        wgt=QWidget()
        wgt.setLayout(layout)
        self.setCentralWidget(wgt)
    def select_dir(self):
        dir_name=QFileDialog.getExistingDirectory(self,"Select Directory","./",QFileDialog.ShowDirsOnly)
        print(dir_name)
        self.dir_name=dir_name
        print(os.listdir(dir_name))
        dirs=glob.glob(dir_name+"/*/")
        print(dirs)
        for dat in dirs:
            print(os.path.basename(dat.rstrip("/")))




if __name__=="__main__":
    app=QApplication(sys.argv)
    win=MainWindow()
    win.show()
    app.exec()
