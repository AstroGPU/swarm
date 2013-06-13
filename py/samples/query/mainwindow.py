from PySide.QtGui import *
from PySide.QtCore import *
from bsddb3.db import DB
from log import LogRecord
from plot import *

class MainWindow(QMainWindow):
    def __init__(self):
        super(MainWindow, self).__init__()
        self.setWindowTitle("Swarm Query")
        openAction = QAction("&Open...", self, shortcut = QKeySequence.Open,
                             triggered=self.openFile)
        exitAction = QAction("E&xit", self, shortcut="Ctrl+Q",
                statusTip="Exit the application",
                triggered=self.close)
        fileMenu = self.menuBar().addMenu("&File")
        fileMenu.addAction(openAction)
        fileMenu.addAction(exitAction)

        self.textEdit = QTextEdit(self)
        self.setCentralWidget(self.textEdit)

    def openFile(self):
        fileName, filtr = QFileDialog.getOpenFileName(self,
             filter="Database files(*.db)")
        if fileName:
            self.loadFile(fileName)

    def loadFile(self,fileName):
        plot_a_file(fileName)
#                s += LogRecord.from_binary(r[1]).bodies_in_keplerian().__str__() + '\n';

#        self.textEdit.setText(s)


