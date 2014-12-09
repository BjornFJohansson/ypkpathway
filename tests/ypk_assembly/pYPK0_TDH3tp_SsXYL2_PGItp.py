#!/usr/bin/env python
# -*- coding: utf-8 -*-
from pydna import read, pcr, parse, Genbank, Assembly2, Dseqrecord

p577,p578,p468,p467,p567,p568 =    parse('''>577        
                                            gttctgatcctcgagcatcttaagaattc                                                
                                            >578          
                                            gttcttgtctcattgccacattcataagt
                                            >468
                                            gtcgaggaacgccaggttgcccact
                                            >467 
                                            ATTTAAatcctgatgcgtttgtctgcacaga
                                            >567
                                            GTcggctgcaggtcactagtgag
                                            >568
                                            GTGCcatctgtgcagacaaacg''')
                                            
from Bio.Restriction import EcoRV

from pYPKpw import pYPKpw

pYPKpw_lin = pYPKpw.linearize(EcoRV)

from pYPKa_Z_TDH3tp import pYPKa_Z_TDH3tp as first
from pYPKa_A_SsXYL2 import pYPKa_A_SsXYL2 as middle
from pYPKa_E_PGItp import pYPKa_E_PGItp as last                                                                               

first  = pcr( p577, p567, first)
middle = pcr( p468, p467, middle)
last   = pcr( p568, p578, last)

asm = Assembly2((pYPKpw_lin, first, middle, last), limit=31)

s = asm.circular_products[0]

pYPK0_TDH3tp_SsXYL2_PGItp = s.synced("tcgcgcgtttcggtgatgacggtgaaaacctctg")

