#!/usr/bin/env python
# -*- coding: utf-8 -*- 
    
from pydna import Dseqrecord, read, pcr

from pYPKa import pYPKa

from Bio.Restriction import ZraI

ins = read("TEF1.txt")

fp = read('''
>pfw579
ttaaatACAATGCATACTTTGTAC
''', ds=False)

rp = read('''
>prv579
taattaaTTTGTAATTAAAACTTAGATTA
''', ds=False)

pYPKa_cut = pYPKa.linearize(ZraI)

ins = pcr(fp, rp, ins)

pYPKa_Z_TEF1tp = (pYPKa_cut + ins).looped().synced(pYPKa)
