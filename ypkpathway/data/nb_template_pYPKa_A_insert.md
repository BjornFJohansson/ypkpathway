#Construction of pYPKa_A_{insert}

Plan for the construction of E. coli vectors [pYPKa_A_{insert}](pYPKa_A_{insert}.gb)
In short, the insert defined below is cloned in pYPKa using the blunt restriction enzyme AjiI

	import pydna

load the vector backbone from a local file.

	pYPKa = pydna.read("pYPKa.gb")

import the restriction enzyme from [Biopython](http://biopython.org/wiki/Main_Page)

	from Bio.Restriction import AjiI

Linearize the vector with both enzymes.

	pYPKa_AjiI  = pYPKa.linearize(AjiI)

Read the insert from a local file.

	ins = pydna.read("{insert}.gb")

Design primers for the terminator promoter. The forward primer adds two adenines to improve
transcription.

	fp_tail = "AA"

	fp, rp = pydna.cloning_primers(ins, fp_tail=fp_tail)

Primers are given the following names:

	fp.id = "{insert}fw"
	rp.id = "{insert}rv"

	print fp.format("fasta")

	print rp.format("fasta")

PCR to create the insert.

	ins = pydna.pcr(fp, rp, ins)

The PCR product has this length in bp.

	len(ins)

Primers annealing on template.

	ins.figure()

A suggested PCR program.

	ins.program()

The final vector is.

	pYPKa_A_{insert} = (pYPKa_AjiI  + ins).looped().synced(pYPKa)

The final vector with reverse insert is.

	pYPKa_A_{insert}b = (pYPKa_AjiI  + ins.rc()).looped().synced(pYPKa)

##Diagnostic PCR confirmation of pYPKa_A_{insert}
The correct structure of each ve is confirmed by 
PCR using primers 577, 342 and the fp defined earlier 
in a multiplex PCR reaction with all three primers present.

	p = {{ x.id: x for x in pydna.parse("primers.fasta") }}

##Expected PCR products sizes:
pYPKa with insert in correct orientation.

    pydna.Anneal( (p['577'], p['342'], fp), pYPKa_A_{insert}).products

pYPKa with insert in reverse orientation.

    pydna.Anneal( (p['577'], p['342'], fp), pYPKa_A_{insert}b).products

Empty pYPKa clone

    pydna.Anneal( (p['577'], p['342'], fp), pYPKa).products

##Diagnostic PCR confirmation 

Calculate cseguid checksum for the resulting plasmids for future reference.

	pYPKa_A_{insert}.cseguid()

The file is named

	pYPKa_A_{insert}.locus = "pYPKa_A_{insert}"[:16]

Stamp sequence with cSEGUID checksum.

	pYPKa_A_{insert}.stamp()

Write sequence to a local file.

	pYPKa_A_{insert}.write("pYPKa_A_{insert}.gb")

#[Download](pYPKa_A_{insert}.gb)

	import pydna
	reloaded = pydna.read("pYPKa_A_{insert}.gb")
	reloaded.verify_stamp()
 
