#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydna import Genbank
gb = Genbank("bjornjobbb@gmail.com")
pCAPs = gb.nucleotide("AJ001614")
pCAPs.features = pCAPs.features [1:]

#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydna import Genbank
gb = Genbank("bjornjobbb@gmail.com")
pSU0 = gb.nucleotide("AB215109")
pSU0.features = pSU0.features[1:]

#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pSU0 import pSU0
from pydna import parse, pcr

primers = parse('''
>470_pSU0f-dEcoRV (50-mer)
gtttactaaaaacacatgtggatattttgactgatttttccatggagggc

>469_pSU0r-dEcoRV (50-mer)
gccctccatggaaaaatcagtcaaaatatccacatgtgtttttagtaaac''', ds=False)

pSU0_EcoRV = pcr( primers, pSU0)[50:].looped().synced("ggatccatcggaattcatattgaaaaagga")

pSU0_EcoRV.add_feature(2741,3545, label= "URA3")


#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pCAPs import pCAPs
from pSU0_EcoRV import pSU0_EcoRV
from pydna import Assembly2
from Bio.Restriction import BamHI, EcoRI, PvuI

psu = pSU0_EcoRV.cut(BamHI, EcoRI).pop(0)

pcaps = pCAPs.cut(PvuI).pop()

a = Assembly2([psu, pcaps], limit=200)

pYPK0 = a.circular_products[0]

pYPK0=pYPK0.synced(pCAPs)

pYPK0.name = "pYKP0"

#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pCAPs import pCAPs

from pydna import pcr, parse


primers = parse('''>568_pCAPsAjiIR (22-mer)
                   GTGCcatctgtgcagacaaacg
                   >567_pCAPsAjiIF (23-mer)
                   GTcggctgcaggtcactagtgag''', ds=False)

pYPKa = pcr(primers, pCAPs).looped().synced(pCAPs)


#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pYPK0 import pYPK0

from pydna import pcr, parse

primers = parse(''' >pYPKpwR
                    GCATGACGTCaccagacgctatgactcacccggacggca

                    >pYPKpwF
                    GCATGATATCttcacaggcggttttcgcacgtacccatg''', ds=False)

seq = pcr(primers, pYPK0)

seq = seq.looped()

pYPKpw = seq.synced(pYPK0)

from Bio.Restriction import EcoRV, FspAI, ZraI

slask = pYPKpw.linearize(ZraI)
slask = pYPKpw.linearize(FspAI)
slask = pYPKpw.linearize(EcoRV)






























'''

{primer_list}

General primers:

577 gttctgatcctcgagcatcttaagaattc
578 gttcttgtctcattgccacattcataagt
468 gtcgaggaacgccaggttgcccact
467 ATTTAAatcctgatgcgtttgtctgcacaga
567 GTcggctgcaggtcactagtgag
568 GTGCcatctgtgcagacaaacg
775 gcggccgctgacTTAAAT
778 ggtaaatccggatTAATTAA
342 CCTTTTTACGGTTCCTGGCCT


                             >-gene-->
            >-TP-->           \     /           >-TP-->                     
             \   /             \   /             \   /
              \ /               \ /               \ /
               |                 |                 |
 577>      775>| 468>       <567 | 568>       <467 |<778     <578    <342
 |||       ||| | |||         ||| | |||         ||| | |||      |||     |||
 --------------Z-----------------A-----------------E---------------------
|              r                 j                 c                     |
|              a                 i                 o                     |
|              I                 I                 R                     |
|                                                  V                     |
|                                                                        |
 --------------------- pYPKa tp gene tp ---------------------------------



'''

