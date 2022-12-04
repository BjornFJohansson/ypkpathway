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
# # Construction of pYPKa_A_{insert}
#
# This notebook describe the construction of the _E. coli_ vector [pYPKa_A_{insert}](pYPKa_A_{insert}.gb).
#
# ![pYPKa_A plasmid](pYPK_A.png "pYPKa_A plasmid")
#
# A part of the [pydna](https://pypi.python.org/pypi/pydna/) package is imported in the code cell below.

# %%
from pydna.readers import read
from pydna.parsers import parse
from pydna.parsers import parse_primers
from pydna.design  import primer_design
from pydna.amplify import pcr
from pydna.amplify import Anneal

# %% [markdown]
# The vector backbone [pYPKa](pYPKa.gb) is read from a local file.

# %%
pYPKa = read("pYPKa.gb")

# %% [markdown]
# The restriction enzyme [AjiI](http://rebase.neb.com/rebase/enz/AjiI.html) is imported from [Biopython](http://biopython.org)

# %%
from Bio.Restriction import AjiI

# %% [markdown]
# The plasmid is linearized with the enzyme.

# %%
pYPKa_AjiI  = pYPKa.linearize(AjiI)

# %% [markdown]
# The insert sequence is read from a [local file]({insert}.gb).

# %%
ins = read("{insert}.gb")

# %% [markdown]
# Primers are needed to PCR amplify the insert. The forward primer adds two adenines in front of the start codon
# which is a feature commonly found in highly expressed _S. cerevisiae_ genes.

# %%
fp_tail = "AA"

# %% [markdown]
# Primers are designed in the code cell below.

# %%
ins = primer_design(ins)
fp = fp_tail + ins.forward_primer
rp = ins.reverse_primer

# %% [markdown]
# The primers are included in the [new_primer.txt](new_primers.txt) list and in the end of the [pathway notebook](pw.ipynb) file.

# %%
print(fp.format("fasta"))
print(rp.format("fasta"))
with open("new_primers.txt", "a+") as f:
    f.write(fp.format("fasta"))
    f.write(rp.format("fasta"))

# %% [markdown]
# The gene is amplifed using the newly designed primers.

# %%
ins = pcr(fp, rp, ins)

# %% [markdown]
# The PCR product has this length in bp.

# %%
len(ins)

# %% [markdown]
# Primers anneals on template in the figure below.

# %%
ins.figure()

# %% [markdown]
# A suggested PCR program.

# %%
ins.program()

# %% [markdown]
# The final vector is:

# %%
pYPKa_A_{insert} = (pYPKa_AjiI  + ins).looped().synced("gttctgatcctcgagcatcttaagaattc")

# %% [markdown]
# The vector with reverse insert is created below. This vector theoretically make up
# fifty percent of the clones. The PCR strategy below is used to identify the correct clones.

# %%
pYPKa_A_{insert}b = (pYPKa_AjiI  + ins.rc()).looped().synced("gttctgatcctcgagcatcttaagaattc")

# %% [markdown]
# A combination of standard primers and the newly designed primers are
# used for the strategy to identify correct clones.
# Standard primers are listed [here](standard_primers.txt).
# The standard primers are read into a dictonary in the code cell below.

# %%
p = {{ x.id: x for x in parse_primers("standard_primers.txt") }}

# %% [markdown]
# ## Diagnostic PCR confirmation of pYPKa_A_{insert}
# The correct structure of pYPKa_A_{insert} is confirmed by PCR using standard primers
# 577 and 342 that are vector specific together with the {insert}fw primer specific for the insert
# in a multiplex PCR reaction with three primers present in the PCR reaction.
#
# Two PCR products are expected if the insert was sucessfully cloned, sizes depending
# on the orientation of the insert.
# If the vector is empty, only one short product is formed.
#
# ## Expected PCR products sizes:
#
# pYPKa_A_{insert} with insert in correct orientation.

# %%
Anneal( (p['577'], p['342'], fp), pYPKa_A_{insert}).products

# %% [markdown]
# pYPKa_A_{insert} with insert in reverse orientation.

# %%
Anneal( (p['577'], p['342'], fp), pYPKa_A_{insert}b).products

# %% [markdown]
# Empty clone

# %%
Anneal( (p['577'], p['342'], fp), pYPKa).products

# %% [markdown]
# The cseguid checksum for the resulting plasmid is calculated for future reference.
# The [cseguid checksum](http://pydna.readthedocs.org/en/latest/pydna.html#pydna.utils.cseguid)
# uniquely identifies a circular double stranded sequence.

# %%
pYPKa_A_{insert}.cseguid()

# %% [markdown]
# The file is given a name based on the cloned insert

# %%
pYPKa_A_{insert}.locus = "pYPKa_A_{insert}"[:16]

# %% [markdown]
# Sequence is stamped with the cseguid checksum.
# This can be used to verify the integrity of the sequence file.

# %%
pYPKa_A_{insert}.stamp("cSEGUID")

# %% [markdown]
# The sequence is written to a local file.

# %%
pYPKa_A_{insert}.write("pYPKa_A_{insert}.gb")
