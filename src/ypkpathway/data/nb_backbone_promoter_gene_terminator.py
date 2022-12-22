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
# # {name}
#
# This notebook describes the assembly of the
# [_Saccharomyces cerevisiae_](http://www.yeastgenome.org)
# transcriptional unit vector.
#
# It is made from a linear vector fragment and three PCR products:
#
# - a linearized {backbone} vector for maintenance in _S. cerevisiae_ or _E. coli_ (red dashed line in figure below)
# - a promoter PCR product from the `{promoter}` vector
# - a gene PCR product from the `{gene}` vector
# - a terminator PCR product from the `{terminator}` vector
#
# The four linear DNA fragments are joined by _in-vivo_ homologous recombination in a
# [_Saccharomyces cerevisiae_](http://wiki.yeastgenome.org/index.php/Commonly_used_strains) strain.
#
# ![tu](tu.png "tu")
#
# A part of the [pydna](https://pypi.python.org/pypi/pydna/) package is imported in the code cell below.

# %%
from pydna.parsers import parse_primers
from pydna.readers import read
from pydna.amplify import pcr
from pydna.assembly import Assembly

# %% [markdown]
# The Yeast Pathway Kit [standard primers](standard_primers.fasta) are read into a dictionary in the code cell below.

# %%
p = {{x.name: x for x in parse_primers("standard_primers.fasta")}}

# %% [markdown]
# The backbone vector [{backbone}]({backbone}.gb) is read from a local file in the code cell below.

# %%
backbone = read("{backbone}.gb")

# %% [markdown]
# The backbone vector is linearized by digestion with [{enz}](https://www.google.com/search?q={enz}).

# %%
from Bio.Restriction import {enz}

# %%
linear_backbone = backbone.linearize({enz})

# %% [markdown]
# The pYPKa derived _E. coli_ plasmids containing 
# - [promoter]({promoter}.gb)
# - [gene]({gene}.gb)
# - [terminator]({terminator}.gb)
#
# are read into three variables below.

# %%
promoter_template   = read("{promoter}.gb")
gene_template       = read("{gene}.gb")
terminator_template = read("{terminator}.gb")

# %% [markdown]
# ### PCR
# Three DNA fragments are PCR amplified using [standard primers](standard_primers.fasta).
#
# [Suggested PCR programs](#Suggested-PCR-programs) can be found at the end of this document.
# %%

fp_prom = p['{fp_prom}']
rp_prom = p['{rp_prom}']
fp_gene = p['{fp_gene}']
rp_gene = p['{rp_gene}']
fp_term = p['{fp_term}']
rp_term = p['{rp_term}']

# %%
prom = pcr(fp_prom, rp_prom, promoter_template)
gene = pcr(fp_gene, rp_gene, gene_template)
term = pcr(fp_term, rp_term, terminator_template)

# %%
prom.name = "{promoter}"[8:]
gene.name = "{gene}"[8:]
term.name = "{terminator}"[8:]

# %% [markdown]
#
# The fragments will be assembled by _in-vivo_ [homologous recombination](http://www.ncbi.nlm.nih.gov/pubmed/2828185):

# %%
asm = Assembly((linear_backbone, prom, gene, term), limit=31)
asm

# %% [markdown]
# The Assembly object above should normally indicate four fragments and eight nodes.

# %%
candidates = asm.assemble_circular()
candidates
# %% [markdown]
# There should normally be two candidates of equal size. These sequences should be identical.

# %%
candidate, *rest = candidates
# %%
candidate.cseguid() == rest[0].cseguid()
# %%
candidate.figure()
# %% [markdown]
# The candidate vector is synchronized to the 577 primer. This means that
# the plasmid origin is shifted so that it matches the backbone vector.

# %%
result = candidate.synced(fp_prom)

# %% [markdown]
# ### Diagnostic PCR confirmation
#
# The structure of the final vector is confirmed by two
# separate PCR reactions, one for the promoter and gene and
# one for the gene and terminator.
#
# PCR using standard primers 577 and 467 to amplify promoter and gene.

# %%
product = pcr( fp_prom, rp_gene, result)

# %% [markdown]
# A correct clone should give this size in base pairs:

# %%
print(len(product))

# %% [markdown]
# If the promoter is missing from the assembly, the PCR product will have this size in base pairs:

# %%
print(len(product) - len(prom))

# %% [markdown]
# If the gene is missing from the assembly, the PCR product will have this size in base pairs:

# %%
print(len(product) - len(gene))

# %% [markdown]
# PCR using standard primers 468 and 578 to amplify gene and terminator.

# %%
product2 = pcr( fp_gene, rp_term, result)

# %% [markdown]
# A correct clone should give this size:

# %%
print(len(product2))

# %% [markdown]
# If the gene is missing from the assembly, the PCR product will have this size in base pairs:

# %%
print(len(product2) - len(gene))

# %% [markdown]
# If the terminator is missing from the assembly, the PCR product will have this size in base pairs:

# %%
print(len(product2) - len(term))

# %% [markdown]
# The cseguid checksum for the resulting plasmid is calculated for future reference.
# The [cseguid checksum](http://pydna.readthedocs.org/en/latest/pydna.html#pydna.utils.cseguid)
# uniquely identifies a circular double stranded sequence.

# %%
print(result.cseguid())

# %% [markdown]
# The file name is based on the promoter, gene and terminator designations.

# %%
result.locus = f"tu:{{gene.name}}"[:16]
result.definition = "{name}"

# %% [markdown]
# Sequence is stamped with cseguid checksum. This can be used to verify the
# integrity of the sequence file.

# %%
result.stamp("cSEGUID")

# %% [markdown]
# Write sequence to a local file.

# %%
result.write(f"{{result.definition}}.gb")

# %% [markdown]
# ### Suggested PCR programs
#
# For the [amplification](#PCR) of promoter, gene and terminator.
# %%
print(prom.forward_primer.name, prom.reverse_primer.name)
print(prom.program())

# %%
print(gene.forward_primer.name, gene.reverse_primer.name)
print(gene.program())

# %%
print(term.forward_primer.name, term.reverse_primer.name)
print(term.program())
