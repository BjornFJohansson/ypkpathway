# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
# ---

# %% [markdown]
# # Pathway {name}
#
# This notebook describes the assembly of {length} single gene expression cassettes into a single pathway.
# Notebooks describing the single gene expression vectors are linked at the end of this document as are notebooks
# describing pYPKa promoter, gene and terminator vectors. Specific primers needed are also listed below.
#
# ![pathway with N genes](pw.png "pathway with N genes")
#
# The assembly performed here is based on content of the [INDATA.txt](INDATA.txt) text file.
# A part of the [pydna](https://pypi.python.org/pypi/pydna/) package is imported in the code cell below.

# %%
from pydna.parsers import parse_primers
from pydna.readers import read
from pydna.amplify import pcr
from pydna.assembly import Assembly

# %% [markdown]
# Initiate the standard primers needed to amplify each cassette.
# The first cassette in the pathway is amplified with standard
# primers 577 and 778, the last with
# 775 and 578 and all others with 775 and 778.
# Standard primers are listed [here](standard_primers.txt).

# %%
p = {{ x.id: x for x in parse_primers("standard_primers.txt") }}

# %% [markdown]
# The backbone vector is linearized with [EcoRV](http://rebase.neb.com/rebase/enz/EcoRV.html).

# %%
from Bio.Restriction import EcoRV, NotI, PacI

pYPKpw = read("pYPKpw.gb")

# %% [markdown]
# The cassette_products variable holds the list of expression cassette PCR products fragments to
# be assembled.

# %%
cassette_products = []

# %% [markdown]
# The expression cassettes comes from a series of single gene expression vectors
# held in the template_vectors list.

# %%
cassette_vectors = ({cas_vectors})
template_vectors = [read(v.strip()) for v in cassette_vectors if v.strip()]
for tv in  template_vectors:
    display(tv)

# %% [markdown]
# The first cassette in the pathway is amplified with standard primers 577 and 778. Suggested PCR conditions can be found at the end of this document.

# %%
cassette_products.append( pcr( p['577'], p['778'],  template_vectors[0] ) )

# %% [markdown]
# Cassettes in the middle cassettes are amplified with standard primers 775 and 778. Suggested PCR conditions can be found at the end of this document.

# %%
cassette_products.extend( pcr( p['775'], p['778'], v) for v in template_vectors[1:-1] )

# %% [markdown]
# The last cassette in the pathway is amplified with standard primers 775 and 578. Suggested PCR conditions can be found at the end of this document.

# %%
cassette_products.append( pcr( p['775'], p['578'], template_vectors[-1] ) )

# %% [markdown]
# The cassettes are given names based on their order in the final construct in the code cell below.

# %%
for i, cp in enumerate(cassette_products):
    cp.name = "Cassette {{}}".format(i+1)
    print(cp.name)

# %% [markdown]
# Cassettes and plasmid backbone are joined by homologous recombination in a Saccharomyces cerevisiae ura3 host
# which selects for the URA3 gene in pYPKpw.

# %%
asm = Assembly( [pYPKpw.linearize(EcoRV)] + cassette_products, limit=167-47-10)
asm

# %% [markdown]
# Normally, only one circular product should be formed since the
# homology limit is quite large (see cell above). More than one
# circular products might indicate an incorrect strategy.
# The largest recombination product is chosen as candidate for
# the {name} pathway.

# %%
candidate = asm.assemble_circular()[0]

# %% [markdown]
# This assembly figure shows how the fragments came together.

# %%
candidate.figure()

# %% [markdown]
# The candidate vector is synchronized to the 577 primer. This means that
# the plasmid origin is shifted so that it matches the backbone vector.

# %%
pw = candidate.synced("gttctgatcctcgagcatcttaagaattc")

# %% [markdown]
# The cseguid checksum for the resulting plasmid is calculated for future reference.
# The [cseguid checksum](http://pydna.readthedocs.org/en/latest/pydna.html#pydna.utils.cseguid)
# uniquely identifies a circular double stranded sequence.

# %%
pw.cseguid()

# %% [markdown]
# The file is given a name based on the sequence of expressed genes.

# %%
pw.locus = "pw"
pw.definition = "{name}"

# %% [markdown]
# Stamp sequence with cseguid checksum. This can be used to verify the
# integrity of the sequence file.

# %%
pw.stamp("cSEGUID")

# %% [markdown]
# Write sequence to a local file.

# %%
pw.write("{name}.gb")

# %% [markdown]
# The pathway can be extended by digestion with either NotI or PacI or both provided that the enzymes cut once in the final pathway sequence.

# %%
print("NotI cuts {{}} time(s) and PacI cuts {{}} time(s) in the final pathway.".format(len(pw.cut(NotI)), len(pw.cut(PacI))))

# %% [markdown]
# ### New Primers needed for assembly.
#
# This list contains all needed primers that are not in the standard primer [list](standard_primers.txt) above.

# %%
try:
    with open("new_primers.txt") as f:
        text = f.read()
except IOError:
    text = "no new primers needed."
print(text)

# %% [markdown]
# ### New single gene expression vectors ( pYPK0_promoter_gene_terminator ) needed for assembly.
#
# Hyperlinks to notebook files describing the single gene expression plasmids needed for the assembly.
#
# {tp_gene_tp_links}
#
# ### New pYPKa vectors needed for assembly of the single gene expression vectors above.
#
# Hyperlinks to notebook files describing the pYPKa plasmids needed for the assembly of the single gene clones listed above.
#
# {pYPKa_clones}
#
# ### Suggested PCR conditions

# %%
for prd in cassette_products:
    print("\n\n")
    print("product name:", prd.name)
    print("forward primer", prd.forward_primer.name)
    print("reverse primer", prd.reverse_primer.name)
    print(prd.program())
