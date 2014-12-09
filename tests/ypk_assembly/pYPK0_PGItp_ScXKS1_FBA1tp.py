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

from pYPKa_Z_PGItp import pYPKa_Z_PGItp as first
from pYPKa_A_ScXKS1 import pYPKa_A_ScXKS1 as middle
from pYPKa_E_FBA1tp import pYPKa_E_FBA1tp as last                                                                               

first  = pcr( p577, p567, first)
middle = pcr( p468, p467, middle)
last   = pcr( p568, p578, last)

asm = Assembly2((pYPKpw_lin, first, middle, last), limit=31)

s = asm.circular_products[0]

pYPK0_PGItp_ScXKS1_FBA1tp = s.synced("tcgcgcgtttcggtgatgacggtgaaaacctctg")

