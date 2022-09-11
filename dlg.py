import sys
from PyQt5.QtWidgets import (QMainWindow, QTextEdit, QAction, QFileDialog, QApplication)
from PyQt5.QtGui import QIcon

class Example(QMainWindow):
    def __init__(self): # default constructor
        super().__init__()
        self.showDialog()
        self.initUI()   # call init User Interface

    def initUI(self):
        self.textEdit=QTextEdit()   # text editor object?
        self.setCentralWidget(self.textEdit) # place 
        input()
        self.statusBar()

        openFile=QAction(QIcon('sample.jpg'),'Open',self)
        openFile.setShortcut('Ctrl+O')
        openFile.setStatusTip("Open new File")

        openFile.triggered.connect(self.showDialog)

        menubar=self.menuBar()
        fileMenu=menubar.addMenu('&File')
        fileMenu.addAction(openFile)

        self.setGeometry(300,300,350,300)
        self.setWindowTitle("File dialog")
        self.show()
        
    def showDialog(self):
        fname=QFileDialog.getOpenFileName(self,"Open file","/home")
        fp=open(fname,"r")
        print(fp.readlines())
        """
        if fname[0]:
            f=open(fname[0],"r")
            with f:
                data=f.read()
                self.textEdit.setText(data)
        """

if __name__=="__main__":
    app=QApplication(sys.argv)
    ex=Example()
    sys.exit(app.exec_())
