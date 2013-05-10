from PySide.QtGui import *
from sys import argv
from mainwindow import MainWindow

a = QApplication(argv)
w = MainWindow()
w.show()
a.exec_()
