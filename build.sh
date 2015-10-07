#!/bin/bash                 # This “shebang” tells what program to use to interpret the script.
conda install bsddb
$PYTHON setup.py build
$PYTHON setup.py install    # Python command to install the script.
