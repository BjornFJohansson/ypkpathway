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

from pYPK0_TEF1tp_SsXYL1_TDH3tp import pYPK0_TEF1tp_SsXYL1_TDH3tp as cas1
from pYPK0_TDH3tp_SsXYL2_PGItp import pYPK0_TDH3tp_SsXYL2_PGItp as cas2
from pYPK0_PGItp_ScXKS1_FBA1tp import pYPK0_PGItp_ScXKS1_FBA1tp as cas3
from pYPK0_FBA1tp_ScTAL1_PDC1tp import pYPK0_FBA1tp_ScTAL1_PDC1tp as cas4
cas1  = pcr( p577, p778, cas1)
cas2  = pcr( p775, p778, cas2)
cas3  = pcr( p775, p778, cas3)
cas4 = pcr( p775, p578, cas4)    

asm = Assembly2((pYPKpw_lin, cas1,cas2,cas3,cas4), limit=167-47-10)

seq = asm.circular_products[0]
        
pYPK0_TEF1tp_SsXYL1_TDH3tp_SsXYL2_PGItp_ScXKS1_FBA1tp_ScTAL1_PDC1tp_pw = seq.synced(pYPKpw)
