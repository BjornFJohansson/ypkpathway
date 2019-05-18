#!/usr/bin/env python
# -*- coding: utf-8 -*-
# editor.py

import syntax
from PyQt5 import  QtGui, QtWidgets



app = QtWidgets.QApplication([])
editor =  QtWidgets.QPlainTextEdit()
highlight = syntax.PythonHighlighter(editor.document())
editor.show()

# Load syntax.py into the editor for demo purposes
infile = open('syntax.py', 'r')
editor.setPlainText(infile.read())

app.exec_()