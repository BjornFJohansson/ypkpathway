#!/usr/bin/env python
# -*- coding: utf-8 -*- 
    
from pydna import Dseqrecord, read, pcr

from pYPKa import pYPKa

from Bio.Restriction import AjiI

ins = read("SsXYL2.txt")

fp = read('''
>pfw1092
aaATGACTGCTAACCCTTC
''', ds=False)

rp = read('''
>prv1092
TTACTCAGGGCCGTCA
''', ds=False)

pYPKa_cut = pYPKa.linearize(AjiI)

ins = pcr(fp, rp, ins)

pYPKa_A_SsXYL2 = (pYPKa_cut + ins).looped().synced(pYPKa)
