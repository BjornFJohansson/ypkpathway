#!/usr/bin/env python
# -*- coding: utf-8 -*-

execfile("header.py")

def design():
    #---------------------------------------------------------------------------
    # This script should define a Dseqrecord named seq
    from pYPK0 import pYPK0

    from pydna import pcr, parse


    primers = parse(''' >pYPKpwR
                        GCATGACGTCaccagacgctatgactcacccggacggca

                        >pYPKpwF
                        GCATGATATCttcacaggcggttttcgcacgtacccatg''', ds=False)

    seq = pcr(primers, pYPK0)

    seq = seq.looped()

    seq = seq.synced(pYPK0)

    from Bio.Restriction import EcoRV, FspAI, ZraI

    slask = seq.linearize(ZraI)
    slask = seq.linearize(FspAI)
    slask = seq.linearize(EcoRV)

    # This script should define a Dseqrecord named seq
    #---------------------------------------------------------------------------
    assert isinstance(seq, Dseqrecord)
    return seq

execfile("footer.py")
