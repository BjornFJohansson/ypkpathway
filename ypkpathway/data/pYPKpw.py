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
