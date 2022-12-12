#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""docstring."""

from pydna.amplify import pcr
from Bio.Restriction import ZraI, EcoRV, AjiI
from pydna.readers import read
from pydna.genbank import genbank
import csv
from pydna.myprimers import PrimerList
p = PrimerList()

pYPKa = read("pYPKa.gb")

pYPKa_Z = pYPKa.linearize(ZraI)
pYPKa_A = pYPKa.linearize(AjiI)
pYPKa_E = pYPKa.linearize(EcoRV)

bbdict = {"Z": pYPKa_Z, "A": pYPKa_A, "E": pYPKa_E}
fdict = {"Z": "promoter", "A": "gene", "E": "terminator"}

csvfile = "promoter_list_003.csv"
with open(csvfile) as csvfile:
    reader = csv.reader(csvfile)
    for row in reader:
        ins, letter, pf, pr, template, comment = row

        temp = genbank(template)

        insert = pcr(p[int(pf)], p[int(pr)], temp)
        
        label = f"{ins}_promoter"
        
        if letter == "Z":
        
            insert.add_feature(type="promoter", 
                               name=label,
                               ApEinfo_label=label,
                               ApEinfo_fwdcolor="#FC6600",
                               ApEinfo_revcolor="#7C4700",
                               fp=p[int(pf)].format("tab").strip(),
                               rp=p[int(pr)].format("tab").strip(),
                               reference=template)

        if letter == "E":
            
            insert.add_feature(type="terminator", 
                               name=label,
                               ApEinfo_label=label,
                               ApEinfo_fwdcolor="#7C4700",
                               ApEinfo_revcolor="#FC6600",
                               fp=p[int(pf)].format("tab").strip(),
                               rp=p[int(pr)].format("tab").strip(),                               
                               reference=template)

        bb = bbdict[letter]

        plasmid = (bb + insert).looped().synced(p[577])
        
        sfs = insert.extract_feature(-1)
        
        for n in range(len(plasmid.features)):
            f = plasmid.extract_feature(n)
            if f.seq == sfs.seq:
                plasmid.features = [plasmid.features[n]] + plasmid.features[:n] + plasmid.features[n+1:]
                break  
                
        pname = {"Z": "pYPKa_Z_", "A": "pYPKa_A_", "E": "pYPKa_E_"}[letter]
        pname += ins
        plasmid.locus = pname
        plasmid.definition = pname
        plasmid.accession = "accession"
        plasmid.annotations["date"] = "13-JUL-2021"
        plasmid.stamp()


        # plasmid.annotations["comment"] += 

        plasmid.write(pname+".gb")

        #from pydna.editor import ape
        #ape(plasmid)        
        #if letter == "E": break
    
        print(pname)
    
    
    