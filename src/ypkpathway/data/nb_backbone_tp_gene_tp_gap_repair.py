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
# # pYPK0_{tpz}_{gene}_{tpe}
#
# This notebook describes the assembly of the [_Saccaromyces cerevisiae_](www.yeastgenome.org)
# single gene expression vector pYPK0_{tpz}_{gene}_{tpe}.
#
# It is made by _in-vivo_ homologous recombination between three PCR products and one linear vector fragment.
# The PCR products are a promoter generated from a pYPK_Z vector, a gene from a template and
# a terminator from a pYPKa_E vector. The three PCR products are joined to
# a linearized [pYPKpw](https://github.com/BjornFJohansson/ypk-xylose-pathways/blob/master/pYPKpw.gb)
# backbone vector that has the [URA3](http://www.yeastgenome.org/locus/S000000747/overview)
# marker and a _S. crevisiae_ [2 micron](http://blog.addgene.org/plasmids-101-yeast-vectors) origin of replication.
#
# The four linear DNA fragments are joined by homologous recombination in a
# [_Saccharomyces cerevisiae_](http://wiki.yeastgenome.org/index.php/Commonly_used_strains) ura3 mutant.
#
# ![pYPK0_promoter_gene_terminator](tp_g_tp.png "pYPK0_promoter_gene_terminator")
#
# A part of the [pydna](https://pypi.python.org/pypi/pydna/) package is imported in the code cell below.

# %%
from pydna.parsers import parse_primers
from pydna.readers import read
from pydna.design import primer_design
from pydna.amplify import pcr
from pydna.assembly import Assembly

# %% [markdown]
# The YPK [standard primers](standard_primers.txt) are read into a dictionary in the code cell below.

# %%
p = {{ x.id: x for x in parse_primers("standard_primers.txt") }}

# %% [markdown]
# The backbone vector [pYPKpw](pYPKpw.gb) is read from a local file in the code cell below.

# %%
pYPKpw = read("pYPKpw.gb")

# %% [markdown]
# The backbone vector is linearized by digestion with [EcoRV](http://rebase.neb.com/rebase/enz/EcoRV.html).
# The restriction enzyme functionality is provided by [biopython](http://biopython.org).

# %%
from Bio.Restriction import EcoRV

pYPK_EcoRV = pYPKpw.linearize(EcoRV)

# %% [markdown]
# The pYPKa derived _E. coli_ plasmids containing the [promoter](pYPKa_Z_{tpz}.gb) and [terminator](pYPKa_E_{tpe}.gb)
# as well as the [gene template]({gene}.gb) sequence are read into three variables in the code cell below.

# %%
promoter_template   = read("pYPKa_Z_{tpz}.gb")
gene_template       = read("{gene}.gb")
terminator_template = read("pYPKa_E_{tpe}.gb")

# %% [markdown]
# The construction of the two vector above are described in the [pYPKa_ZE_{tpz}](pYPKa_ZE_{tpz}.ipynb) notebooks.
#
# The promoter is amplified with from [pYPKa_Z_{tpz}](pYPKa_Z_{tpz}.gb). A suggested PCR program can be found at the end of this document.

# %%
prom = pcr( p['577'], p['567'], promoter_template)

# %% [markdown]
# Primers with tails are needed for the recombination between the gene and the promoter/terminator PCR products.
# The tails are designed to provide 33 bp of terminal homology between the DNA fragments.

# %%
fp_tail = "tgcccactttctcactagtgacctgcagccgacAA"
rp_tail = "AAatcctgatgcgtttgtctgcacagatggCAC"

# %% [markdown]
# Primers with the tails above are designed for the gene template in the code cell below.

# %%
ins = primer_design(gene_template)
fp = fp_tail + ins.forward_primer
rp = rp_tail + ins.reverse_primer

# %% [markdown]
# The primers are included in the [new_primers.txt](new_primers.txt) list and in the end of the [pathway notebook](pw.ipynb) file.

# %%
print(fp.format("fasta"))
print(rp.format("fasta"))
with open("new_primers.txt", "a+") as f:
    f.write(fp.format("fasta"))
    f.write(rp.format("fasta"))

# %% [markdown]
# The gene is amplifed using the newly designed primers. A suggested PCR program can be found at the end of this document.

# %%
gene = pcr( fp, rp, gene_template)

# %% [markdown]
# The terminator is amplified from [pYPKa_E_{tpe}](pYPKa_E_{tpe}.gb). A suggested PCR program can be found at the end of this document.

# %%
term = pcr( p['568'], p['578'], terminator_template)

# %% [markdown]
# The four linear DNA fragments are mixed and transformed
# to a _Saccharomyces cerevisiae_ ura3 mutant.
#
# The fragments will be assembled by _in-vivo_ [homologous recombination](http://www.ncbi.nlm.nih.gov/pubmed/2828185):

# %%
asm = Assembly( (pYPK_EcoRV, prom, gene, term), limit=31 )

asm

# %% [markdown]
# The representation of the asm object above should normally indicate one circcular product only.
# More than one circular products might indicate an incorrect assembly strategy or represent
# by-products that might arise in the assembly process.
# The largest recombination product is chosen as candidate for the pYPK0_{tpz}_{gene}_{tpe} vector.

# %%
candidate = asm.assemble_circular()[0]

candidate.figure()

# %% [markdown]
# The candidate vector is synchronized to the 577 primer. This means that
# the plasmid origin is shifted so that it matches the backbone vector.

# %%
result = candidate.synced("gttctgatcctcgagcatcttaagaattc")

# %% [markdown]
# ### Diagnostic PCR confirmation
#
# The structure of the final vector is confirmed by two
# separate PCR reactions, one for the promoter and gene and
# one for the gene and terminator.
#
# PCR using standard primers 577 and 467 to amplify promoter and gene.

# %%
product = pcr( p['577'], p['467'], result)

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
product2 = pcr( p['468'], p['578'], result)

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
result.cseguid()

# %% [markdown]
# The file is named based on the names of promoter, gene and terminator.

# %%
result.locus = "pYPK0_tp_g_tp"
result.definition = "pYPK0_{tpz}_{gene}_{tpe}"

# %% [markdown]
# Sequence is stamped with cseguid checksum. This can be used to verify the
# integrity of the sequence file.

# %%
result.stamp("cSEGUID")

# %% [markdown]
# Write sequence to a local file.

# %%
result.write("pYPK0_{tpz}_{gene}_{tpe}.gb")

# %% [markdown]
# ## PCR programs for the amplification of Promoter, Gene and Terminator
#
# see cell # 6
#
# Promoter

# %%
prom.program()

# %% [markdown]
# Gene

# %%
gene.program()

# %% [markdown]
# Terminator

# %%
term.program()
