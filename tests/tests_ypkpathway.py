#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil
import nose
import textwrap
import pydna
import tempfile
from ypkpathway import pathway

tmp = os.path.join( tempfile.gettempdir(), "ypkpathway_test_dir" )

def test_ypk():

    datafiles = '''pth1.txt|pYPK0_CiGXF1_PsXYL2.gb|iM8oDuvJPMPO995IdW3B0oo0Hkc
                   pth2.txt|pYPK0_SsXYL1_SsXYL2.gb|CB_qLhPgemW0XNLQOQEAdJKFujU
                   pth3.txt|pYPK0_NC_006038_CiGXF1_PsXYL2.gb|48_Bek9U1wxXlq1otmE7YHjYpnk
                   pth4.txt|pYPK0_SsXYL1_SsXYL2_ScXKS1.gb|5OxynmwQA3br0cKAG8It7VVNGrg
                   pth5.txt|pYPK0_SsXYL1_SsXYL2_ScXKS1_ScTAL1.gb|K8z4ijkYa0hA0KEOhv7-6PNJgBM
                   pth6.txt|pYPK0_SsXYL1_SsXYL2_ScXKS1_ScTAL1.gb|K8z4ijkYa0hA0KEOhv7-6PNJgBM
                   pth7.txt|pYPK0_SsXYL1_SsXYL2_ScXKS1_ScTAL1.gb|K8z4ijkYa0hA0KEOhv7-6PNJgBM'''

    for datafile in textwrap.dedent(datafiles).split():

        try:
            shutil.rmtree(tmp)
        except OSError:
            pass

        file_, name, code = datafile.split("|")
        print
        print "############################"
        print file_
        print "############################"

        with open(file_, "rU",) as f: text = f.read()

        pw = pathway( pydna.parse(text), tmp)

        s = pydna.read( os.path.join(tmp, name) )

        assert( s.seguid() == code )

    try:
        shutil.rmtree(tmp)
    except OSError:
        pass

if __name__ == '__main__':
    nose.runmodule()

