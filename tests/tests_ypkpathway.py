#!/usr/bin/env python
# -*- coding: utf-8 -*-

import nose
import textwrap

import pydna

from ypkpathway import PathWay

def test_ypk():
    datafiles = '''./pth1.txt|pYPK0_RPL12Btp_CiGXF1_TDH3tp_PsXYL2_PGI1tp_pw.txt|iM8oDuvJPMPO995IdW3B0oo0Hkc
                   ./pth2.txt|pYPK0_TEF1tp_SsXYL1_TDH3tp_SsXYL2_PGItp_pw.txt|CB/qLhPgemW0XNLQOQEAdJKFujU
                   ./pth3.txt|pYPK0_RPL12Atp_NC_006038_RPL12Btp_CiGXF1_TDH3tp_PsXYL2_PGI1tp_pw.txt|48/Bek9U1wxXlq1otmE7YHjYpnk
                   ./pth4.txt|pYPK0_TEF1tp_SsXYL1_TDH3tp_SsXYL2_PGItp_ScXKS1_FBA1tp_pw.txt|5OxynmwQA3br0cKAG8It7VVNGrg
                   ./pth5.txt|pYPK0_TEF1tp_SsXYL1_TDH3tp_SsXYL2_PGItp_ScXKS1_FBA1tp_ScTAL1_PDC1tp_pw.txt|K8z4ijkYa0hA0KEOhv7+6PNJgBM
                   ./pth6.txt|pYPK0_TEF1tp_SsXYL1_TDH3tp_SsXYL2_PGItp_ScXKS1_FBA1tp_ScTAL1_PDC1tp_pw.txt|K8z4ijkYa0hA0KEOhv7+6PNJgBM'''

    for datafile in textwrap.dedent(datafiles).split():

        file_, name, code = datafile.split("|")

        print file_

        with open(file_, "rU",) as f: text = f.read()

        pw = PathWay( text )
        pw.generate_files()
        s = pydna.read(pw.files[name])
        assert( s.seguid() == code )

if __name__ == '__main__':
    nose.runmodule()

