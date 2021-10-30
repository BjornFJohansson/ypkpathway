#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import shutil
import pytest
import textwrap

from pydna.readers import read

import tempfile
from ypkpathway import pathway

tmp = os.path.join( tempfile.gettempdir(), "ypkpathway_test_dir" )

def test_ypk():

    datafiles = '''pth1.txt|pYPK0_RPL12B_CiGXF1_TDH3_PsXYL2_PGI1.gb|iM8oDuvJPMPO995IdW3B0oo0Hkc
                   pth2.txt|pYPK0_TEF1_SsXYL1_TDH3_SsXYL2_PGI.gb|CB_qLhPgemW0XNLQOQEAdJKFujU
                   pth3.txt|pYPK0_RPL12A_NC_006038_RPL12B_CiGXF1_TDH3_PsXYL2_PGI1.gb|48_Bek9U1wxXlq1otmE7YHjYpnk
                   pth4.txt|pYPK0_TEF1_SsXYL1_TDH3_SsXYL2_PGI_ScXKS1_FBA1.gb|5OxynmwQA3br0cKAG8It7VVNGrg
                   pth5.txt|pYPK0_TEF1_SsXYL1_TDH3_SsXYL2_PGI_ScXKS1_FBA1_ScTAL1_PDC1.gb|K8z4ijkYa0hA0KEOhv7-6PNJgBM
                   pth6.txt|pYPK0_TEF1_SsXYL1_TDH3_SsXYL2_PGI_ScXKS1_FBA1_ScTAL1_PDC1.gb|K8z4ijkYa0hA0KEOhv7-6PNJgBM
                   pth7.txt|pYPK0_TEF1_SsXYL1_TDH3_SsXYL2_PGI_ScXKS1_FBA1_ScTAL1_PDC1.gb|K8z4ijkYa0hA0KEOhv7-6PNJgBM'''

    for pYPKa_A in (True, False):
        
        for datafile in textwrap.dedent(datafiles).split():

            file_, name, code = datafile.split("|")

            print()
            print("############################")
            print("datafile = ", file_)
            print("pYPKa_A  = ", pYPKa_A)
            print("############################")

            with open(file_, "r",) as f: text = f.read()

            try:
                shutil.rmtree(tmp)
            except OSError:
                pass

            pathway( text, tmp, pYPKa_A=pYPKa_A)

            s = read( os.path.join(tmp, name) )

            with open(code+".txt") as f: c = f.read()

            assert "".join( x for x in c.lower() if not x.isspace()) == str(s.seq).lower()



if __name__ == '__main__':
    pytest.cmdline.main([__file__, "-v", "-s"])

