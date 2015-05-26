#Construction of vectors pYPKa_Z_{tp} and pYPKa_E_{tp}

Plan for the construction of E. coli vectors [pYPKa_Z_{tp}](pYPKa_Z_{tp}.gb) and [pYPKa_E_{tp}](pYPKa_E_{tp}.gb)

In short, the insert defined below is cloned in pYPKa using the blunt restriction enzymes ZraI or EcoRV

	import pydna

load the vector backbone from a local file.

	pYPKa = pydna.read("pYPKa.gb")

import the restriction enzymes from [Biopython](http://biopython.org/wiki/Main_Page)

	from Bio.Restriction import ZraI, EcoRV

Linearize the vector with both enzymes.

	pYPKa_ZraI  = pYPKa.linearize(ZraI)
	pYPKa_EcoRV = pYPKa.linearize(EcoRV)

Read the insert from a local file.

	ins = pydna.read("{tp}.gb")

Design primers for the terminator promoter. The primers has specific tails.

	fp_tail = "ttaaat"
	rp_tail = "taattaa"

	fp, rp = pydna.cloning_primers(ins, fp_tail=fp_tail, rp_tail=rp_tail)

Primers are given the following names:

	fp.id = "{tp}fw"
	rp.id = "{tp}rv"

	print fp.format("fasta")

	print rp.format("fasta")

PCR to create the insert.

	prd = pydna.pcr(fp, rp, ins)

The PCR product has this length in bp.

	len(prd)

Primers annealing on template.

	prd.figure()

A suggested PCR program.

	prd.program()

The final vectors are

	pYPKa_Z_{tp} = (pYPKa_ZraI  + prd).looped().synced(pYPKa)
	pYPKa_E_{tp} = (pYPKa_EcoRV + prd).looped().synced(pYPKa)

The final vectors with reverse inserts are.

	pYPKa_Z_{tp}b = (pYPKa_ZraI  + prd.rc()).looped().synced(pYPKa)
	pYPKa_E_{tp}b = (pYPKa_EcoRV + prd.rc()).looped().synced(pYPKa)

##Diagnostic PCR confirmation of pYPKa_Z_{tp}

The correct structure of pYPKa_Z_{tp} is confirmed by PCR using primers
577, 342 and the forward primer that was constructed in a multiplex PCR
 reaction with all three primers present.

	p = {{ x.id: x for x in pydna.parse("primers.fasta") }}

####Expected PCR products sizes:

pYPKa with insert in correct orientation.

    pydna.Anneal( (p['577'], p['342'], fp), pYPKa_Z_{tp}).products

pYPKa with insert in reverse orientation.

    pydna.Anneal( (p['577'], p['342'], fp), pYPKa_Z_{tp}b).products

Empty pYPKa clone

    pydna.Anneal( (p['577'], p['342'], fp), pYPKa).products

##Diagnostic PCR confirmation of pYPKa_E_{tp}

The correct structure of pYPKa_Z_{tp} is confirmed by PCR using primers
568, 342 and the forward primer that was constructed in a multiplex PCR
 reaction with all three primers present.

####Expected PCR products sizes:

pYPKa with insert in correct orientation.

    pydna.Anneal( (p['577'], p['342'], fp), pYPKa_E_{tp}).products

pYPKa with insert in reverse orientation.

    pydna.Anneal( (p['577'], p['342'], fp), pYPKa_Z_{tp}b).products

Empty pYPKa clone

    pydna.Anneal( (p['577'], p['342'], fp), pYPKa).products

Calculate cseguid checksum for the resulting plasmids for future reference.

	print pYPKa_Z_{tp}.cseguid()
	print pYPKa_E_{tp}.cseguid()

The sequences are named.

	pYPKa_Z_{tp}.locus = "pYPKa_Z_{tp}"[:16]
	pYPKa_E_{tp}.locus = "pYPKa_Z_{tp}"[:16]

Stamp sequence with cSEGUID checksum.

	pYPKa_Z_{tp}.stamp()
	pYPKa_E_{tp}.stamp()

Write sequence to a local file.

	pYPKa_Z_{tp}.write("pYPKa_Z_{tp}.gb")
	pYPKa_E_{tp}.write("pYPKa_E_{tp}.gb")

#[Download](pYPKa_Z_{tp}.gb)

	import pydna
	reloaded = pydna.read("pYPKa_Z_{tp}.gb")
	reloaded.verify_stamp()

#[Download](pYPKa_Z_{tp}.gb)

	import pydna
	reloaded = pydna.read("pYPKa_E_{tp}.gb")
	reloaded.verify_stamp()

