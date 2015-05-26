#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import nose

os.environ["PYTHONDONTWRITEBYTECODE"] = "True"
sys.dont_write_bytecode = True
os.environ["pydna_cache"]  = "nocache"

def main():
    cwd = os.getcwd()
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(os.path.join(dname,"tests"))
    print "Python version:", sys.version
    print "Operating system:", os.name, sys.platform
    print "running tests"
    print
    nose.run(argv=[__file__, "--all-modules", "--verbosity=3", "--nocapture", "--with-doctest", "--doctest-options=+ELLIPSIS"])
    #os.chdir(os.path.join(dname, "ypkpathway"))
    #print "doctests"
    #nose.run(argv=[__file__, "--all-modules", "--verbosity=3", "--nocapture", "--with-doctest", "--doctest-options=+ELLIPSIS"])
    os.chdir(cwd)

if __name__ == '__main__':
    main()