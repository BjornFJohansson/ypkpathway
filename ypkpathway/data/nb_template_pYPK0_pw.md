#Pathway {name}

Import the [pydna](https://pypi.python.org/pypi/pydna/) functionality.

    import pydna

Initiate the primers needed to amplify each cassette.
The first cassette in the pathway is amplified with 
primers 577 and 778, the last with
775 and 578 and all others with 775 and 778.

    p = {{ x.id: x for x in pydna.parse("primers.fasta") }}

    from Bio.Restriction import EcoRV

    pYPKpw = pydna.read("pYPKpw.gb")

The backbone vector is linearized with EcoRV.

    assembly_fragments = [ pYPKpw.linearize(EcoRV) ]

The expression cassettes comes from a series of single gene expression vectors.

    cas_vectors ='''
{cas_vectors}'''.splitlines()

    template_vectors = [pydna.read(v.strip()) for v in cas_vectors if v.strip()]

    template_vectors

The first cassette in the pathway is amplified with primers 577 and 778

    assembly_fragments.append( pydna.pcr( p['577'], p['778'],  template_vectors[0] ) )

Middle cassettes are amplified with 775 and 778.

    assembly_fragments.extend( pydna.pcr( p['775'], p['778'], v) for v in template_vectors[1:-1] ) 

The last cassette in the pathway is amplified with primers 775 and 578

    assembly_fragments.append( pydna.pcr( p['775'], p['578'], template_vectors[-1] ) )

Cassettes and plasmid backbone are joined by homologous recombination in a Saccharomyces cerevisiae ura3 host.

    asm = pydna.Assembly( assembly_fragments, limit=167-47-10)
    asm

The largest recombination product is chosen.

    candidate = asm.circular_products[0]

Assembly figure.
            
    candidate.figure()

The final pathway is synced to the backbone vector.

    pw = candidate.synced(pYPKpw)

Calculate cseguid checksum for the resulting plasmid for future reference.

    pw.cseguid()

The file is named.

    pw.locus = "pw"
    pw.definition = "{name}"

Stamp sequence with cSEGUID checksum.

    pw.stamp()

Write sequence to a local file.

    pw.write("{name}.gb")

#[Download]({name}.gb)

    import pydna

    reloaded = pydna.read("{name}.gb")

    reloaded.verify_stamp()

###New Primers needed for assembly.

{primer_list}

###New pYPKa clones needed for assembly.

{pYPKa_clones}

###New pYPK0_tp_gene_tp clones needed for assembly.

{tp_gene_tp_links}
	

