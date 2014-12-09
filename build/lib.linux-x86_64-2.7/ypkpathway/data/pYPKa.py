#!/usr/bin/env python
# -*- coding: utf-8 -*- 

execfile("header.py")

def design():
    #---------------------------------------------------------------------------
    # This script should define a Dseqrecord named seq
    
    from pCAPs import pCAPs
    
    from pydna import pcr, parse
    
    
    primers = parse('''>568_pCAPsAjiIR (22-mer)
                       GTGCcatctgtgcagacaaacg
                       >567_pCAPsAjiIF (23-mer)
                       GTcggctgcaggtcactagtgag''', ds=False)
    
    seq = pcr(primers, pCAPs)

    seq = seq.looped()
    
    seq = seq.synced(pCAPs)    
    
    # This script should define a Dseqrecord named seq
    #---------------------------------------------------------------------------
    assert isinstance(seq, Dseqrecord)
    return seq
    
execfile("footer.py")
