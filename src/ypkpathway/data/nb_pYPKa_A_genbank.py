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
from pydna.genbank import Genbank
from pydna.parsers import parse_primers
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
# Access to [Genbank](http://www.ncbi.nlm.nih.gov/nuccore) is needed in order to download the template.
# If the email address below is not yours, change it before executing this script as you must always give NCBI a way to contact you when using their service.

# %%
gb = Genbank("{email}")

# %% [markdown]
# The template is downloaded from Genbank below.

# %%
template  = gb.nucleotide("{gblink}")
template

# %% [markdown]
# The two primers below are used to amplify the insert.

# %%
fp,rp =  parse_primers(""">{fpn}
                          {fps}
                          >{rpn}
                          {rps}""")

# %% [markdown]
# The gene is amplifed using the primers specified above.

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