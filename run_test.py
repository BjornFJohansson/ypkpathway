#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import pytest

os.environ["PYTHONDONTWRITEBYTECODE"] = "True"
sys.dont_write_bytecode = True
os.environ["pydna_cache"]  = "nocache"

def main():
    cwd = os.getcwd()
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(os.path.join(dname,"tests"))
    print("Python version:", sys.version)
    print("Operating system:", os.name, sys.platform)
    print("running tests")
    print()
    
    try:
        import coveralls
    except ImportError:
        print("python-coveralls NOT installed!")
        args = []
    else:
        del coveralls
        args = ["--cov=pydna", "--cov-report=html", "--cov-report=xml"]    
    try:
        import nbval
    except ImportError:
        print("nbval NOT installed!")
    else:
        del nbval
        args.append("--nbval")   

    args = [".", "-v", "-s"] + args 
    cwd = os.getcwd()
    pytest.cmdline.main(args)
    os.chdir(cwd)
    
    
    os.chdir(cwd)

if __name__ == '__main__':
    main()