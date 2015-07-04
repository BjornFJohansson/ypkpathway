import sys
from PyQt4 import QtGui #, QtCore
#from PyQt4.QtCore import Qt

from pkg_resources import resource_filename

class Main(QtGui.QMainWindow):

    def __init__(self, parent = None):
        QtGui.QMainWindow.__init__(self,parent)

        self.filename = ""

        self.initUI()

    def assemble(self):
        print self.text.toPlainText()

    def initToolbar(self):

        self.openAction = QtGui.QAction(QtGui.QIcon( resource_filename("ypkpathway","icons/open.png")),"Open file",self)
        self.openAction.setStatusTip("Open existing document")
        self.openAction.setShortcut("Ctrl+O")
        self.openAction.triggered.connect(self.open)

        self.saveAction = QtGui.QAction(QtGui.QIcon( resource_filename("ypkpathway","icons/save.png")),"Save",self)
        self.saveAction.setStatusTip("Save document")
        self.saveAction.setShortcut("Ctrl+S")
        self.saveAction.triggered.connect(self.save)

        self.cutAction = QtGui.QAction(QtGui.QIcon( resource_filename("ypkpathway","icons/cut.png")),"Cut to clipboard",self)
        self.cutAction.setStatusTip("Delete and copy text to clipboard")
        self.cutAction.setShortcut("Ctrl+X")
        self.cutAction.triggered.connect(self.text.cut)

        self.copyAction = QtGui.QAction(QtGui.QIcon( resource_filename("ypkpathway","icons/copy.png")),"Copy to clipboard",self)
        self.copyAction.setStatusTip("Copy text to clipboard")
        self.copyAction.setShortcut("Ctrl+C")
        self.copyAction.triggered.connect(self.text.copy)

        self.pasteAction = QtGui.QAction(QtGui.QIcon( resource_filename("ypkpathway","icons/paste.png")),"Paste from clipboard",self)
        self.pasteAction.setStatusTip("Paste text from clipboard")
        self.pasteAction.setShortcut("Ctrl+V")
        self.pasteAction.triggered.connect(self.text.paste)

        self.undoAction = QtGui.QAction(QtGui.QIcon( resource_filename("ypkpathway","icons/undo.png")),"Undo last action",self)
        self.undoAction.setStatusTip("Undo last action")
        self.undoAction.setShortcut("Ctrl+Z")
        self.undoAction.triggered.connect(self.text.undo)

        self.redoAction = QtGui.QAction(QtGui.QIcon( resource_filename("ypkpathway","icons/redo.png")),"Redo last undone thing",self)
        self.redoAction.setStatusTip("Redo last undone thing")
        self.redoAction.setShortcut("Ctrl+Y")
        self.redoAction.triggered.connect(self.text.redo)

        self.assembleAction = QtGui.QAction(QtGui.QIcon( resource_filename("ypkpathway","icons/assemble.png")),"Start Assembly",self)
        self.assembleAction.setStatusTip("Start Assembly")
        self.assembleAction.triggered.connect(self.assemble)

        self.toolbar = self.addToolBar("Options")

        self.toolbar.addAction(self.openAction)
        self.toolbar.addAction(self.saveAction)

        self.toolbar.addSeparator()

        self.toolbar.addAction(self.cutAction)
        self.toolbar.addAction(self.copyAction)
        self.toolbar.addAction(self.pasteAction)
        self.toolbar.addAction(self.undoAction)
        self.toolbar.addAction(self.redoAction)

        self.toolbar.addSeparator()
        self.toolbar.addAction(self.assembleAction)


        # Makes the next toolbar appear underneath this one
        self.addToolBarBreak()



    def initFormatbar(self):

      self.formatbar = self.addToolBar("Format")


    def initMenubar(self):

      menubar = self.menuBar()

      file = menubar.addMenu("File")
      edit = menubar.addMenu("Edit")
      assemble = menubar.addMenu("Start Assembly")

      file.addAction(self.openAction)
      file.addAction(self.saveAction)

      edit.addAction(self.undoAction)
      edit.addAction(self.redoAction)
      edit.addAction(self.cutAction)
      edit.addAction(self.copyAction)
      edit.addAction(self.pasteAction)

      assemble.addAction(self.assembleAction)

    def initUI(self):

        self.text = QtGui.QTextEdit(self)

        self.text.setFontFamily("inconsolata")
        self.text.setFontPointSize(12.0)

        self.initToolbar()
        self.initFormatbar()
        self.initMenubar()

        # Set the tab stop width to around 33 pixels which is
        # about 8 spaces
        self.text.setTabStopWidth(33)

        self.setCentralWidget(self.text)

        # Initialize a statusbar for the window
        self.statusbar = self.statusBar()

        # If the cursor position changes, call the function that displays
        # the line and column number
        self.text.cursorPositionChanged.connect(self.cursorPosition)

        # x and y coordinates on the screen, width, height
        self.setGeometry(0,0,1500,800)

        self.setWindowTitle("ypkpathway GUI")

        self.setWindowIcon(QtGui.QIcon( resource_filename("ypkpathway","icons/ypkpathway.png")))

    def open(self):

        # Get filename and show only .writer files
        self.filename = QtGui.QFileDialog.getOpenFileName(self, 'Open File',".")

        if self.filename:
            with open(self.filename,"rt") as file:
                self.text.setText(file.read())

    def save(self):

        # Only open dialog if there is no filename yet
        if not self.filename:
          self.filename = QtGui.QFileDialog.getSaveFileName(self, 'Save File')

        # Append extension if not there yet
        if not self.filename.endswith(".writer"):
          self.filename += ".writer"

        # We just store the contents of the text file along with the
        # format in html, which Qt does in a very nice way for us
        with open(self.filename,"wt") as file:
            file.write(self.text.toHtml())


    def cursorPosition(self):

        cursor = self.text.textCursor()

        # Mortals like 1-indexed things
        line = cursor.blockNumber() + 1
        col = cursor.columnNumber()

        self.statusbar.showMessage("Line: {} | Column: {}".format(line,col))


def main():

    app = QtGui.QApplication(sys.argv)

    main = Main()
    main.show()

    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
