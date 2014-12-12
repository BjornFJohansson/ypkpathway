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
