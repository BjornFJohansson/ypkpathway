#!/usr/bin/env python
# -*- coding: utf-8 -*- 
from pydna import read, pcr, parse, Assembly2, Dseqrecord

(p577, p578, p775, p778) = parse('''  >577          
                                      gttctgatcctcgagcatcttaagaattc                                                                           
                                      >578            
                                      gttcttgtctcattgccacattcataagt
                                      >775 
                                      gcggccgctgacTTAAAT
                                      >778 
                                      ggtaaatccggatTAATTAA''', ds=False)

from Bio.Restriction import EcoRV

from pYPKpw import pYPKpw

pYPKpw_lin = pYPKpw.linearize(EcoRV)

cas1 = read('pYPK0_RPL12Btp_CiGXF1_TDH3tp.txt')
cas2 = read('pYPK0_TDH3tp_PsXYL2_PGI1tp.txt')
cas1  = pcr( p577, p778, cas1)
cas2 = pcr( p775, p578, cas2)    

asm = Assembly2((pYPKpw_lin, cas1,cas2), limit=167-47-10)

seq = asm.circular_products[0]
        
pYPK0_RPL12Btp_CiGXF1_TDH3tp_PsXYL2_PGI1tp_pw = seq.synced(pYPKpw)
