#!/usr/bin/env python
# -*- coding: utf-8 -*- 
    
from pydna import Dseqrecord, read, pcr

from pYPKa import pYPKa

from Bio.Restriction import EcoRV

ins = read("PDC1.txt")

fp = read('''
>pfw955
ttaaatAGGGTAGCCTCCCCAT
''', ds=False)

rp = read('''
>prv955
taattaaTTTGATTGATTTGACTGT
''', ds=False)

pYPKa_cut = pYPKa.linearize(EcoRV)

ins = pcr(fp, rp, ins)

pYPKa_E_PDC1tp = (pYPKa_cut + ins).looped().synced(pYPKa)
