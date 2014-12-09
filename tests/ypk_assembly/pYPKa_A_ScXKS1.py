#!/usr/bin/env python
# -*- coding: utf-8 -*- 
    
from pydna import Dseqrecord, read, pcr

from pYPKa import pYPKa

from Bio.Restriction import AjiI

ins = read("ScXKS1.txt")

fp = read('''
>pfw1803
aaATGTTGTGTTCAGTAATTCA
''', ds=False)

rp = read('''
>prv1803
TTAGATGAGAGTCTTTTCCA
''', ds=False)

pYPKa_cut = pYPKa.linearize(AjiI)

ins = pcr(fp, rp, ins)

pYPKa_A_ScXKS1 = (pYPKa_cut + ins).looped().synced(pYPKa)
