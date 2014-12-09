#!/usr/bin/env python
# -*- coding: utf-8 -*- 

execfile("header.py")

def design():
    #---------------------------------------------------------------------------
    # This script should define a Dseqrecord named seq

    from pCAPs import pCAPs
    from pSU0_EcoRV import pSU0_EcoRV
    from pydna import pcr, parse, Genbank, Assembly
    from Bio.Restriction import ZraI, AjiI, EcoRV, BamHI, EcoRI, PvuI
    
    psu = pSU0_EcoRV.cut(BamHI, EcoRI).pop(0)
    
    pcaps = pCAPs.cut(PvuI).pop()
    
    a = Assembly([psu, pcaps])
    
    print a.analyze_overlaps(limit=200)
    
    print a.create_graph()
    
    print a.assemble_hr_circular()

    seq = a.circular_products[0]
    
    seq=seq.synced(pCAPs)
    
    seq.name = "pYKP0"

    
    # This script should define a Dseqrecord named seq
    #---------------------------------------------------------------------------
    assert isinstance(seq, Dseqrecord)
    return seq
    
execfile("footer.py")
