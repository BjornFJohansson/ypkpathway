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
# # Construction of pYPK0_A_{insert}
#
# This notebook describe the construction of the _E. coli_ vector [pYPK0_A_{insert}](pYPK0_A_{insert}.gb).
# Primers needed for the amplification of the insert are designed in this notebook.
#
# ![pYPK0_A plasmid](pYPK0_A.png "pYPK0_A plasmid")
#
# A part of the [pydna](https://pypi.python.org/pypi/pydna/) package is imported in the code cell below.

# %%
from pydna.readers import read
from pydna.genbank import Genbank
from pydna.parsers import parse_primers
from pydna.design  import primer_design
from pydna.amplify import pcr
from pydna.assembly import Assembly
from pydna.amplify import Anneal

# %% [markdown]
# The vector backbone [pYPK0](pYPK0.gb) is read from a local file.

# %%
pYPK0 = read("pYPK0.gb")

# %% [markdown]
# The restriction enzyme [AscI](http://rebase.neb.com/rebase/enz/AscI.html) is imported from [Biopython](http://biopython.org)

# %%
from Bio.Restriction import AscI

# %% [markdown]
# The plasmid is linearized with the enzyme.

# %%
pYPK0_AscI  = pYPK0.linearize(AscI)

# %% [markdown]
# Access to [Genbank](http://www.ncbi.nlm.nih.gov/nuccore) is needed in order to download the template.
# If the email address below is not yours, change it before executing this script as you must always give
# NCBI a way to contact you when using their service.

# %%
gb = Genbank("{email}")

# %% [markdown]
# The template is downloaded from Genbank below.

# %%
template  = gb.nucleotide("{gblink}")
template

# %% [markdown]
# Primers are needed to PCR amplify the insert. The forward primer adds two adenines in front of the start codon
# which is a feature commonly found in highly expressed _S. cerevisiae_ genes.

# %%
fp_tail = "actttctcactagtgacctgcagccgacAATG"
rp_tail = "atcctgatgcgtttgtctgcacagatggCAC"

# %% [markdown]
# Primers are designed in the code cell below.

# %%
ins = primer_design(template)
fp = fp_tail + ins.forward_primer
rp = rp_tail + ins.reverse_primer

# %% [markdown]
# The primers are included in the [new_primer.txt](new_primers.txt) file.

# %%
print(fp.format("fasta"))
print(rp.format("fasta"))
with open("new_primers.txt", "a+") as f:
    f.write(fp.format("fasta"))
    f.write(rp.format("fasta"))

# %% [markdown]
# The newly designed primers are used to amplify the insert.

# %%
ins = pcr(fp, rp, template)

# %% [markdown]
# The primers anneal on the template like this.

# %%
ins.figure()

# %% [markdown]
# A suggested PCR program.

# %%
ins.program()

# %% [markdown]
# The linearzed vector and the insert are joined by homologous recombination.

# %%
asm = Assembly((pYPK0_AscI,ins))
asm

# %% [markdown]
# Usually two equally sized products are formed.

# %%
circular_products = asm.assemble_circular()
circular_products

# %% [markdown]
# The first sequence is chosen.

# %%
candidate = circular_products[0]

# %% [markdown]
# The final vector is:

# %%
pYPK0_A_{insert} = candidate.synced("gttctgatcctcgagcatcttaagaattc")

# %% [markdown]
# A combination of standard primers and the gene specific primers are
# used for the strategy to identify correct clones.
# Standard primers are listed [here](standard_primers.txt).
# The standard primers are read into a dictonary in the code cell below.

# %%
p = {{ x.id: x for x in parse_primers("standard_primers.txt") }}

# %% [markdown]
# ## Diagnostic PCR confirmation of pYPK0_A_{insert}
#
# The correct structure of pYPK0_A_{insert} is confirmed by PCR using standard primers
# 577 and 342 that are vector specific together with the {insert}fw primer specific for the insert
# in a multiplex PCR reaction with three primers present in the PCR reaction.
#
# Two PCR products are expected if the insert was sucessfully cloned, sizes depending
# on the orientation of the insert.
# If the vector is empty, only one short product is formed.
#
# ## Expected PCR products sizes:
#
# pYPK0_A_{insert} with insert in correct orientation.

# %%
Anneal( (p['577'], p['342'], fp), pYPK0_A_{insert}).products

# %% [markdown]
# Empty clone

# %%
Anneal( (p['577'], p['342'], fp), pYPK0).products

# %% [markdown]
# The cseguid checksum for the resulting plasmid is calculated for future reference.
# The [cseguid checksum](http://pydna.readthedocs.org/en/latest/pydna.html#pydna.utils.cseguid)
# uniquely identifies a circular double stranded sequence.

# %%
pYPK0_A_{insert}.cseguid()

# %% [markdown]
# The file is given a name based on the cloned insert

# %%
pYPK0_A_{insert}.locus = "pYPK0_A_{insert}"[:16]

# %% [markdown]
# Sequence is stamped with the cseguid checksum.
# This can be used to verify the integrity of the sequence file.

# %%
pYPK0_A_{insert}.stamp("cSEGUID")

# %% [markdown]
# The sequence is written to a local file.

# %%
pYPK0_A_{insert}.write("pYPK0_A_{insert}.gb")
