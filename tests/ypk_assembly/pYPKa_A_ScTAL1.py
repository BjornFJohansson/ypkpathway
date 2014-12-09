#!/usr/bin/env python
# -*- coding: utf-8 -*- 
    
from pydna import Dseqrecord, read, pcr

from pYPKa import pYPKa

from Bio.Restriction import AjiI

ins = read("ScTAL1.txt")

fp = read('''
>pfw1008
aaATGTCTGAACCAGCTC
''', ds=False)

rp = read('''
>prv1008
TTAAGCGGTAACTTTCTT
''', ds=False)

pYPKa_cut = pYPKa.linearize(AjiI)

ins = pcr(fp, rp, ins)

pYPKa_A_ScTAL1 = (pYPKa_cut + ins).looped().synced(pYPKa)
