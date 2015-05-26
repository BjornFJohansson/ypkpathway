# pYPK0_{tpz}_{gene}_{tpe}

This notebook describes the construction of the E. coli/S. cerevisiae
expression vector pYPK0_{tpz}_{gene}_{tpe}.

It is made by _in-vivo_ homologous recombination between
four linear DNA fragments in a Saccharomyces ura3 mutant. 

Import the pydna functionality.

    import pydna

Establish the needed [primers](primers.fasta)

    p = {{ x.id: x for x in pydna.parse("primers.fasta") }}

The backbone vector is [pYPKpw](pYPKpw.gb)

    pYPKpw = pydna.read("pYPKpw.gb")

The backbone vector will be digested by [EcoRV](http://rebase.neb.com/rebase/enz/EcoRV.html)

    from Bio.Restriction import EcoRV

    pYPK_EcoRV = pYPKpw.linearize(EcoRV)

    promoter_template   = pydna.read("pYPKa_Z_{tpz}.gb")
    gene_template       = pydna.read("pYPKa_A_{gene}.gb")
    terminator_template = pydna.read("pYPKa_E_{tpe}.gb")

    prom = pydna.pcr( p['577'], p['567'], promoter_template)
    gene = pydna.pcr( p['468'], p['467'], gene_template)
    term = pydna.pcr( p['568'], p['578'], terminator_template)

The four linear DNA fragments are mixed and transformed
to a Saccharomyces cerevisiae ura3 mutant
The fragments will be assembled by in-vivo homologous recombination:

    asm = pydna.Assembly( (pYPK_EcoRV, prom, gene, term), limit=31 )

    asm

The largest recombination product is chosen.

    candidate = asm.circular_products[0]

    candidate.figure()

The final vector is synchronized with the backbone vector.

    result = candidate.synced(pYPKpw)

###Diagnostic PCR confirmation

The structure of the final vector is confirmed by two
separate PCR reactions, one for the promoter and gene and
one for the gene and terminator.

PCR using primers primer 577 and 467 to amplify promoter and gene.

    product = pydna.pcr( p['577'], p['467'], result)

    print len(product)

    print len(product) - len(prom)

    print len(product) - len(gene)

PCR using primers primer 468 and 578 to amplify gene and terminator.

    product2 = pydna.pcr( p['468'], p['578'], result)

    print len(product2)

    print len(product2) - len(gene)

    print len(product2) - len(term)

Calculate cseguid checksum for the resulting plasmid for future reference.

    result.cseguid()

The file is named.

	result.locus = "pYPK0_tp_g_tp"
    result.definition = "pYPK0_{tpz}_{gene}_{tpe}"

Stamp sequence with cSEGUID checksum.

    result.stamp()

Write sequence to a local file.

    result.write("pYPK0_{tpz}_{gene}_{tpe}.gb")

#[Download](pYPK0_{tpz}_{gene}_{tpe}.gb)

    import pydna
    reloaded = pydna.read("pYPK0_{tpz}_{gene}_{tpe}.gb")
    reloaded.verify_stamp()


