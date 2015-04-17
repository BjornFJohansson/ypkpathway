#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pCAPs import pCAPs

from pydna import pcr, parse


primers = parse('''>568_pCAPsAjiIR (22-mer)
                   GTGCcatctgtgcagacaaacg
                   >567_pCAPsAjiIF (23-mer)
                   GTcggctgcaggtcactagtgag''', ds=False)

pYPKa = pcr(primers, pCAPs).looped().synced(pCAPs)

