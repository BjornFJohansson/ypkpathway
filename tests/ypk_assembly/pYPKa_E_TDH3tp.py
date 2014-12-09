#!/usr/bin/env python
# -*- coding: utf-8 -*- 
    
from pydna import Dseqrecord, read, pcr

from pYPKa import pYPKa

from Bio.Restriction import EcoRV

ins = read("TDH3.txt")

fp = read('''
>pfw698
ttaaatATAAAAAACACGCTTTTTC
''', ds=False)

rp = read('''
>prv698
taattaaTTTGTTTGTTTATGTGTGTT
''', ds=False)

pYPKa_cut = pYPKa.linearize(EcoRV)

ins = pcr(fp, rp, ins)

pYPKa_E_TDH3tp = (pYPKa_cut + ins).looped().synced(pYPKa)
