#!/usr/bin/env python
# -*- coding: utf-8 -*- 

execfile("header.py")

def design():
    #---------------------------------------------------------------------------
    # This script should define a Dseqrecord named seq
    from pSU0 import pSU0
    from pydna import parse, pcr
    
    primers = parse('''
>470_pSU0f-dEcoRV (50-mer)
gtttactaaaaacacatgtggatattttgactgatttttccatggagggc

>469_pSU0r-dEcoRV (50-mer)
gccctccatggaaaaatcagtcaaaatatccacatgtgtttttagtaaac''', ds=False)

    seq = pcr( primers, pSU0)[50:].looped().synced("ggatccatcggaattcatattgaaaaagga")
    
    # This script should define a Dseqrecord named seq
    #---------------------------------------------------------------------------
    assert isinstance(seq, Dseqrecord)
    return seq
    
execfile("footer.py")
