#!/usr/bin/env python
# -*- coding: utf-8 -*- 
    
from pydna import Dseqrecord, read, pcr

from pYPKa import pYPKa

from Bio.Restriction import AjiI

ins = read("SsXYL1.txt")

fp = read('''
>pfw957
aaATGCCTTCTATTAAGTTGAA
''', ds=False)

rp = read('''
>prv957
TTAGACGAAGATAGGAATCTT
''', ds=False)

pYPKa_cut = pYPKa.linearize(AjiI)

ins = pcr(fp, rp, ins)

pYPKa_A_SsXYL1 = (pYPKa_cut + ins).looped().synced(pYPKa)
