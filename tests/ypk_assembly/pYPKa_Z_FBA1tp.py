#!/usr/bin/env python
# -*- coding: utf-8 -*- 
    
from pydna import Dseqrecord, read, pcr

from pYPKa import pYPKa

from Bio.Restriction import ZraI

ins = read("FBA1.txt")

fp = read('''
>pfw630
ttaaatATAACAATACTGACAGTACTAAA
''', ds=False)

rp = read('''
>prv630
taattaaTTTGAATATGTATTACTTGGT
''', ds=False)

pYPKa_cut = pYPKa.linearize(ZraI)

ins = pcr(fp, rp, ins)

pYPKa_Z_FBA1tp = (pYPKa_cut + ins).looped().synced(pYPKa)
