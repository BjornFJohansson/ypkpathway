#Construction of pYPKa_A_{insert}

This notebook describe the construction of the E. coli vector [pYPKa_A_{insert}](pYPKa_A_{insert}.gb)
with an insert for which PCR primers are also designed.

The insert defined below is cloned in pYPKa using the blunt restriction 
enzyme [AjiI](http://rebase.neb.com/rebase/enz/AjiI.html). The insert can be 
cloned between a promoter and a terminator in a single gene expression construct.

![pYPKa_A plasmid](pYPK_A.png "pYPKa_A plasmid")

The [pydna](https://pypi.python.org/pypi/pydna/) package is imported in the code cell below. 
There is a [publication](http://www.biomedcentral.com/1471-2105/16/142) describing pydna as well as
[documentation](http://pydna.readthedocs.org/en/latest/) available online. 
Pydna is developed on [Github](https://github.com/BjornFJohansson/pydna). 

	import pydna

The vector backbone pYPKa is read from a local [file](pYPKa.gb).

	pYPKa = pydna.read("pYPKa.gb")

import the restriction enzyme from [Biopython](http://biopython.org/wiki/Main_Page)

	from Bio.Restriction import AjiI

Linearize the vector.

	pYPKa_AjiI  = pYPKa.linearize(AjiI)

The insert sequence is read from a local file. This sequence was parsed from the ypkpathway data file.

	ins = pydna.read("{insert}.gb")

Primers are designed for the insert. The forward primer adds two adenines in front of the start codon to improve
transcription.

	fp_tail = "AA"

The pydna cloning_primers function is used to design the primers.

	fp, rp = pydna.cloning_primers(ins, fp_tail=fp_tail)

Primers are given the names below. These primers are included in the primer list in the end of the [pathway notebook](pw.ipynb) file.

	fp.id = "{insert}fw"
	rp.id = "{insert}rv"

	print(fp.format("tab"))

	print(rp.format("tab"))

PCR to create the insert.

	ins = pydna.pcr(fp, rp, ins)

The PCR product has this length in bp.

	len(ins)

Primers annealing on template.

	ins.figure()

A suggested PCR program.

	ins.program()

The final vector is:

	pYPKa_A_{insert} = (pYPKa_AjiI  + ins).looped().synced(pYPKa)

The final vector with reverse insert is created below. This vector theoretically make up
fifty percent of the clones. The PCR strategy below is used to identify the correct clones.

	pYPKa_A_{insert}b = (pYPKa_AjiI  + ins.rc()).looped().synced(pYPKa)

A combination of standard primers and the newly designed primers are 
used for the strategy to identify correct clones.
Standard primers are listed [here](primers.fasta).

	p = {{ x.id: x for x in pydna.parse("primers.fasta") }}

##Diagnostic PCR confirmation of pYPKa_A_{insert}
The correct structure of pYPKa_A_{insert} is confirmed by PCR using standard primers
577 and 342 that are vector specific together with the {insert}fw primer specific for the insert 
in a multiplex PCR reaction with all three primers present.

Two PCR products are expected if the insert was cloned, the sizes depend 
on the orientation. If the vector is empty or contains another insert, only one
product is formed.

##Expected PCR products sizes:

pYPKa_A_{insert} with insert in correct orientation.

    pydna.Anneal( (p['577'], p['342'], fp), pYPKa_A_{insert}).products

pYPKa_A_{insert} with insert in reverse orientation.

    pydna.Anneal( (p['577'], p['342'], fp), pYPKa_A_{insert}b).products

Empty clone

    pydna.Anneal( (p['577'], p['342'], fp), pYPKa).products

Calculate cseguid checksum for the resulting plasmids for future reference.

	pYPKa_A_{insert}.cseguid()

The file is given a name based on the cloned insert

	pYPKa_A_{insert}.locus = "pYPKa_A_{insert}"[:16]

Stamp sequence with cseguid checksum. This can be used to verify the 
integrity of the sequence file.

	pYPKa_A_{insert}.stamp()

The sequence is written to a local file.

	pYPKa_A_{insert}.write("pYPKa_A_{insert}.gb")

#[Download](pYPKa_A_{insert}.gb)

	import pydna
	reloaded = pydna.read("pYPKa_A_{insert}.gb")
	reloaded.verify_stamp()
 
