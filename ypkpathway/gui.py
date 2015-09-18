#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import print_function
import sys
import os
import platform
import subprocess
import shutil
import codecs

from PyQt4.QtCore import QRect, QMetaObject, QObject, QEvent
from PyQt4.QtGui  import (QApplication, QMainWindow, QWidget,
                           QGridLayout, QTabWidget, QPlainTextEdit, QCheckBox,
                           QMenuBar, QMenu, QStatusBar, QAction,
                           QIcon, QFileDialog, QMessageBox, QFont)

import qrc
import pydna
import ypkpathway

from pkg_resources import resource_filename




class Assembler(QMainWindow):
    def __init__(self, parent=None):
        super(Assembler, self).__init__(parent)
        self.resize(800, 600)
        self.filename  = None
        self.filetuple = None
        self.dirty = False  # Refers to Data Page only.
        self.nb = None
        centralwidget = QWidget(self)
        gridLayout = QGridLayout(centralwidget)
        self.tabWidget = QTabWidget(centralwidget)


        # textbox

        self.tab = QWidget()
        font = QFont()
        font.setFamily("Inconsolata")
        font.setPointSize(14)
        self.tab.setFont(font)
        gridLayout_3 = QGridLayout(self.tab)
        self.plainTextEdit = QPlainTextEdit(self.tab)
        self.plainTextEdit.installEventFilter(self)
        self.plainTextEdit.setAcceptDrops(True)
        gridLayout_3.addWidget(self.plainTextEdit, 0, 0, 1, 1)
        self.tabWidget.addTab(self.tab, "")





        self.tab_2 = QWidget()
        self.tab_2.setFont(font)
        gridLayout_2 = QGridLayout(self.tab_2)
        self.plainTextEdit_2 = QPlainTextEdit(self.tab_2)
        gridLayout_2.addWidget(self.plainTextEdit_2, 0, 0, 1, 1)
        self.tabWidget.addTab(self.tab_2, "")

        self.tab_3 = QWidget()
        self.tab_3.setFont(font)
        gridLayout_3 = QGridLayout(self.tab_3)
        self.checkbox = QCheckBox("Cloning genes by tailed primers (no pYPKa_A vectors constructed)")
        self.checkbox.setChecked(True)
        gridLayout_3.addWidget(self.checkbox, 0, 0, 0, 0)
        self.tabWidget.addTab(self.tab_3, "")

        gridLayout.addWidget(self.tabWidget, 0, 0, 1, 1)
        self.setCentralWidget(centralwidget)
        menubar = QMenuBar(self)
        menubar.setGeometry(QRect(0, 0, 800, 29))
        menu_File = QMenu(menubar)
        self.menu_Solve = QMenu(menubar)
        self.menu_Help = QMenu(menubar)
        self.setMenuBar(menubar)
        self.statusbar = QStatusBar(self)
        self.setStatusBar(self.statusbar)
        self.action_New = QAction(self)
        self.actionSave_As = QAction(self)
        self.action_Save = QAction(self)
        self.action_Open = QAction(self)
        self.action_Quit = QAction(self)
        self.action_About = QAction(self)
        self.actionShow_CCPL = QAction(self)
        self.action_Solve  = QAction(self)
        self.action_OpenNB = QAction(self)
        self.action_CCPL = QAction(self)
        self.action_Help = QAction(self)
        menu_File.addAction(self.action_New)
        menu_File.addAction(self.action_Open)
        menu_File.addAction(self.actionSave_As)
        menu_File.addAction(self.action_Save)
        menu_File.addSeparator()
        menu_File.addAction(self.action_Quit)
        self.menu_Solve.addAction(self.action_Solve)
        self.menu_Solve.addAction(self.action_OpenNB)
        self.menu_Help.addAction(self.action_About)
        #self.menu_Help.addAction(self.action_CCPL)
        #self.menu_Help.addAction(self.action_Help)
        menubar.addAction(menu_File.menuAction())
        menubar.addAction(self.menu_Solve.menuAction())
        menubar.addAction(self.menu_Help.menuAction())

        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab),\
                                   "Data Page")
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_2),\
                                   "Assembly log")
        self.tabWidget.setTabText(self.tabWidget.indexOf(self.tab_3),\
                                   "Settings")
        menu_File.setTitle("&File")
        self.menu_Solve.setTitle("&Assemble")
        self.menu_Help.setTitle("&About")
        self.tabWidget.setCurrentIndex(0)
        self.action_New.setText("&New")
        self.action_Open.setText("&Open")
        self.actionSave_As.setText("Save &As")
        self.action_Save.setText("&Save")
        self.action_Quit.setText("&Quit")
        self.action_Solve.setText("&Assemble")
        self.action_OpenNB.setText("&Open &pathway")
        self.action_About.setText("&About")
        #self.action_CCPL.setText("&CCPL")
        #self.action_Help.setText("&Help")
        self.action_Quit.triggered.connect(self.close)
        allToolBar = self.addToolBar("AllToolBar")
        allToolBar.setObjectName("AllToolBar")
        self.addActions(allToolBar, (self.action_Open,
                                     self.actionSave_As,
                                     self.action_Save,
                                     self.action_Solve,
                                     self.action_OpenNB,
                                     self.action_Quit ))
        self.action_New.triggered.connect(self.fileNew)
        self.action_Open.triggered.connect(self.fileOpen)
        self.actionSave_As.triggered.connect(self.fileSaveAs)
        self.action_Save.triggered.connect(self.fileSave)
        self.action_Solve.triggered.connect(self.solveAssembly)
        self.action_OpenNB.triggered.connect(self.openNB)
        self.action_About.triggered.connect(self.aboutBox)
        #self.action_CCPL.triggered.connect(self.displayCCPL)
        #self.action_Help.triggered.connect(self.help)
        self.plainTextEdit.textChanged.connect(self.setDirty)
        self.action_New = self.editAction(self.action_New, None,\
                            'ctrl+N', 'filenew', 'New File.')
        self.action_Open = self.editAction(self.action_Open, None,
                            'ctrl+O', 'fileopen', 'Open File.')
        self.actionSave_As = self.editAction(self.actionSave_As,\
                            None, 'ctrl+A', 'filesaveas',\
                            'Save and Name File.')
        self.action_Save = self.editAction(self.action_Save, None,
                            'ctrl+S', 'filesave', 'Save File.')
        self.action_Solve = self.editAction(self.action_Solve, None,
                            '', 'solve', 'Assemble.')
        self.action_OpenNB = self.editAction(self.action_OpenNB, None,
                            '', 'ipynb', 'Open pathway.')
        self.action_About = self.editAction(self.action_About, None,
                            'ctrl+B', 'about','Pop About Box.')
        self.action_CCPL = self.editAction(self.action_CCPL, None,
                            'ctrl+G', 'licence', 'Show Licence')
        self.action_Help = self.editAction(self.action_Help, None,
                            'ctrl+H', 'help', 'Show Help Page.')
        self.action_Quit =  self.editAction(self.action_Quit, None,
                            'ctrl+Q', 'quit', 'Quit the program.')
        self.plainTextEdit_2.setReadOnly(True)

        self.setWindowTitle("ypkpathway")
        self.setWindowIcon(QIcon( resource_filename("ypkpathway","icons/ypkpathway.png")))
        self.plainTextEdit.setFocus()

    def eventFilter(self, object, event):
        #print(event.type(), QEvent.DragEnter, object, self.plainTextEdit)
        if (object is self.plainTextEdit):
            if (event.type() == QEvent.DragEnter):
                if event.mimeData().hasUrls():
                    event.accept()   # must accept the dragEnterEvent or else the dropEvent can't occur !!!
                    print("accept")
                else:
                    event.ignore()
                    print("ignore")
            if (event.type() == QEvent.Drop):
                if event.mimeData().hasUrls():   # if file or link is dropped
                    urlcount = len(event.mimeData().urls())  # count number of drops
                    url = event.mimeData().urls()[0]   # get first url
                    object.setPlainText('abc')   # assign first url to editline
                    event.accept()  # doesnt appear to be needed
                    print(456)
                    return True
            return False # lets the event continue to the edit
        return False

    def setDirty(self):
        '''On change of text in textEdit window, set the flag
        "dirty" to True'''
        index = self.tabWidget.currentIndex()
        if index is not 0:
            return
        if self.dirty:
            return True
        self.dirty = True
        self.updateStatus('self.dirty set to True')

    def clearDirty(self):
        'Clear dirty flag'
        self.dirty = False

    def fileNew(self):
        '''Clear both Data Page and Solution Page.'''
        self.plainTextEdit.setPlainText(' ')
        self.plainTextEdit_2.setPlainText(' ')
        self.clearDirty()
        self.filename = None

    def okToContinue(self):
        if self.dirty:
            reply = QMessageBox.question(self,
                    "Data Loader - Unsaved Changes",
                    "Save unsaved changes?",
                    QMessageBox.Yes|QMessageBox.No|QMessageBox.Cancel)
            if reply == QMessageBox.Cancel:
                return False
            elif reply == QMessageBox.Yes:
                self.clearDirty()
                return self.fileSave()
        return True

    def okRead(self):
        'Pop-up a warning message.'
        reply = QMessageBox.warning(self,
                "Warning",
                '''\nFile Open and Save only in Data Page
\n\(Use SaveAs for the Assembly log)''', QMessageBox.Ok)
        return True

    def fileOpen(self):
        '''Open a file in Data Page (with index == 0)'''
        if self.tabWidget.currentIndex():
            self.okRead()
            return
        if not self.okToContinue():
            return
        dir_ = (os.path.dirname(unicode(self.filename)) if self.filename is not None else ".")
        self.filename = unicode(QFileDialog.getOpenFileName(self,"Open File", dir_,))
        #print(unicode(self.filename))
        #self.filename = self.filetuple[0]
        fname = self.filename
        #  QFileDialog returns a tuple x with x[0] = file name and
        #  x[1] = type of filter.
        if fname:
            self.loadFile(fname)
            self.updateStatus('New file opened.')

    def loadFile(self, fname=None):
        fl = codecs.open(fname, "U", "utf8")
        text = fl.read()
        self.plainTextEdit.setPlainText(text)
        self.dirty = False

    def fileSave(self):
        '''Save file with current file name.'''
        if self.tabWidget.currentIndex():
            self.okRead()
            return
        if self.filename is None:
            return self.fileSaveAs()
        else:
            flname = self.filename
            if flname:
                tempText = self.plainTextEdit.toPlainText()
                with codecs.open(flname, 'w', 'utf8') as fl: fl.write(tempText)
                self.dirty = False
                self.updateStatus('File saved.')
                return True
            else:
                self.updateStatus('Failed to save... ')
                return False
        self.filename = None
        self.dirty = False

    def fileSaveAs(self):
        '''Save file with a new name.'''
        qpr = self.qprintline
        fname = self.filename or "NoName.txt"
        self.filename = unicode(QFileDialog.getSaveFileName(self,"ypkpathway - Save File", fname))
        flname = self.filename or "NoName.txt"
        self.filename = flname
        fl = codecs.open(flname, 'wU', 'utf8')
        tempText = unicode(self.plainTextEdit.toPlainText())
        fl.write(tempText)
        fl.close()
        self.dirty = False
        self.updateStatus('File saved.')


    def solveAssembly(self):
        printline = self.qprintline
        self.plainTextEdit_2.clear()
        self.tabWidget.setCurrentIndex(1)
        flbase = os.path.basename(unicode(self.filename))

        title = 'Assembly log for ' + flbase

        printline('='*len(title))
        printline(title)
        printline('='*len(title))


        qstringobj = self.plainTextEdit.toPlainText()
        #print(type(qstringobj)) #<class 'PyQt4.QtCore.QString'>
        #print(qstringobj.toUtf8()[3268:3279])
        #print(str(qstringobj.toUtf8()[3268:3279]))
        #print(type(rawtext), "rawtext")
        #codec0 = .QTextCodec.codecForName("UTF-16");
        #rawtext = unicode(codec0.fromUnicode(tmp), 'UTF-16')
        #unicode(qstringobj.toUtf8(), encoding="UTF-8").decode()

        pth = pydna.parse(str(qstringobj.toUtf8()))

        if len(pth)==0:
            printline("No of sequences found in Data window")
            return

        if self.filename is None:
            self.fileSaveAs()

        dir_, ext = os.path.splitext( unicode(self.filename))

        fl, log = ypkpathway.pathway( pth, dir_, pYPKa_A = not self.checkbox.isChecked(), print = printline)

        if not fl:
            return

        with codecs.open(os.path.join(dir_, u"log.txt"),"wU","utf8") as f: f.write(log)

        shutil.copy2( unicode(self.filename), os.path.join(dir_, u"INDATA_"+os.path.basename(unicode(self.filename))))

        printline('')
        printline('\n\nAssembly finished.')
        printline('click on the Open pathway button above to open the pathway in the default web browser')
        self.nb =  fl.path

    def qprintline(self, line):
        '''Append one line to Solution Page.'''
        self.plainTextEdit_2.appendPlainText(line.rstrip()) #.decode("utf8"))
        QApplication.processEvents()

    def openNB(self):
        if self.nb:
            subprocess.Popen(["ipython", "notebook", self.nb])



    def aboutBox(self):

        from PyQt4.QtCore import QT_VERSION_STR
        from PyQt4.Qt import PYQT_VERSION_STR
        from sip import SIP_VERSION_STR
        from _version import get_versions
        __version__ = get_versions()["version"][:5]
        del get_versions
        from IPython import __version__ as IPython_version

        QMessageBox.about(self, "About ypkassembler",
                             u"""<b>Planning yeast pathway kit constructions.</b>
                                 <p>Copyright 2015 Björn Johansson.
                                 This software is released under a BSD style license.
                                 This software comes with no warranties
                                 expressed or implied.<br><br>
                                 Python version: {}<br><br>
                                 ypkassembler version: {}<br>
                                 IPython version: {}<br>
                                 Qt version: {}<br>
                                 SIP version: {}<br>
                                 PyQt version: {}<br>
                                 pydna version: {}<br>
                                 """.format(sys.version,
                                            __version__,
                                            IPython_version,
                                            QT_VERSION_STR,
                                            SIP_VERSION_STR,
                                            PYQT_VERSION_STR,
                                            pydna.__version__[:5]))




    def displayCCPL(self):
        '''Read and display CCPL licence.'''
        self.plainTextEdit.setPlainText(open('CCPL.txt').read())
        self.dirty = False
        self.filename = 'COPYING.txt'
        self.updateStatus('CCPL displayed.')

    def help(self):
        '''Read and display a help file- currently the README.txt.'''
        self.plainTextEdit.setPlainText(open('README.md').read())
        self.dirty = False
        self.filename = 'README.txt'
        self.updateStatus('README displayed.')

    def addActions(self, target, actions):
        '''Actions are added to Tool Bar.'''
        for action in actions:
            if action is None:
                target.addSeparator()
            else:
                target.addAction(action)

    def editAction(self, action, slot=None, shortcut=None, icon=None,
                     tip=None):
        '''This method adds to action: icon, shortcut, ToolTip,\
        StatusTip and can connect triggered action to slot '''
        if icon is not None:
            action.setIcon(QIcon(":/%s.png" % (icon)))
        if shortcut is not None:
            action.setShortcut(shortcut)
        if tip is not None:
            action.setToolTip(tip)
            action.setStatusTip(tip)
        if slot is not None:
            action.triggered.connect(slot)
        return action

    def qreadline(self, lineNo):
        '''Read one line from Data Page (lineNo starts with 0)'''
        return unicode(self.plainTextEdit.document().\
            findBlockByLineNumber(lineNo).text()).rstrip()



    def updateStatus(self, message):
        '''Keep status current.'''
        if self.filename is not None:
            flbase = os.path.basename(unicode(self.filename))
            self.setWindowTitle(unicode("ypkpathway - " +\
                                         flbase + "[*]") )
            self.statusBar().showMessage(message, 5000)
            self.setWindowModified(self.dirty)


def main():
    app = QApplication(sys.argv)
    frame = Assembler()
    frame.show()
    sys.exit(app.exec_())

if __name__ == "__main__":
    main()
