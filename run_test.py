#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import pytest

def main():
    cwd = os.getcwd()
    os.chdir("tests")
    args = [".", "-v", "-s"]
    pytest.cmdline.main(args)
    os.chdir(cwd)
    
if __name__ == '__main__':
    main()