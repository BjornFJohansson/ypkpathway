#!/usr/bin/env python
# -*- coding: utf-8 -*- 
    
from pydna import Dseqrecord, read, pcr

from pYPKa import pYPKa

from Bio.Restriction import EcoRV

ins = read("PGI.txt")

fp = read('''
>pfw1000
ttaaatAATTCAGTTTTCTGACTGA
''', ds=False)

rp = read('''
>prv1000
taattaaTTTTAGGCTGGTATCTTG
''', ds=False)

pYPKa_cut = pYPKa.linearize(EcoRV)

ins = pcr(fp, rp, ins)

pYPKa_E_PGItp = (pYPKa_cut + ins).looped().synced(pYPKa)
