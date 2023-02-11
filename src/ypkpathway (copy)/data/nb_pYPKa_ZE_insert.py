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
# # Construction of promoter vector pYPKa_Z_{tp} and terminator vector pYPKa_E_{tp}
#
# This notebook describe the construction of _E. coli_ vectors [pYPKa_Z_{tp}](pYPKa_Z_{tp}.gb) and [pYPKa_E_{tp}](pYPKa_E_{tp}.gb)
# with the same insert for which PCR primers are also designed.
#
# The insert defined below is cloned in pYPKa using the blunt restriction
# enzymes [ZraI](http://rebase.neb.com/rebase/enz/ZraI.html) and [EcoRV](http://rebase.neb.com/rebase/enz/EcoRV.html) in
# two different plasmids. The insert cloned in [ZraI](http://rebase.neb.com/rebase/enz/ZraI.html)
# will be used as a promoter, while in the [EcoRV](http://rebase.neb.com/rebase/enz/EcoRV.html) site the insert will be used as a
# terminator.
#
# ![pYPKa_Z and pYPKa_E plasmids](pYPK_ZE.png "pYPKa_Z and pYPKa_E plasmids")
#
# A part of the [pydna](https://pypi.python.org/pypi/pydna/) package is imported in the code cell below.

# %%
from pydna.readers import read
from pydna.parsers import parse
from pydna.parsers import parse_primers
from pydna.design import primer_design
from pydna.amplify import pcr
from pydna.amplify import Anneal

# %% [markdown]
# The vector backbone [pYPKa](pYPKa.gb) is read from a local file.

# %%
pYPKa = read("pYPKa.gb")

# %% [markdown]
# Both restriction enzymes are imported from [Biopython](http://biopython.org)

# %%
from Bio.Restriction import ZraI, EcoRV

# %% [markdown]
# The vector is linearized with both enzymes.

# %%
pYPKa_ZraI  = pYPKa.linearize(ZraI)
pYPKa_EcoRV = pYPKa.linearize(EcoRV)

# %% [markdown]
# The insert sequence is read from a local file. This sequence was parsed from the ypkpathway data file.

# %%
ins = read("{tp}.gb")

# %% [markdown]
# Primers for the terminator promoter need specific tails in order to produce
# a [SmiI](http://rebase.neb.com/rebase/enz/SmiI.html) and a [PacI](http://rebase.neb.com/rebase/enz/PacI.html)
# when cloned in pYPKa in the EcoRV cloning position.

# %%
fp_tail = "ttaaat"
rp_tail = "taattaa"

# %% [markdown]
# Primers with the tails above are designed in the code cell below.

# %%
ins = primer_design(ins)
fp = fp_tail + ins.forward_primer
rp = rp_tail + ins.reverse_primer

# %% [markdown]
# The primers are included in the [new_primer.txt](new_primers.txt) list and in the end of the [pathway notebook](pw.ipynb) file.

# %%
print(fp.format("fasta"))
print(rp.format("fasta"))
with open("new_primers.txt", "a+") as f:
    f.write(fp.format("fasta"))
    f.write(rp.format("fasta"))

# %% [markdown]
# PCR to create the insert using the newly designed primers.

# %%
prd = pcr(fp, rp, ins)

# %% [markdown]
# The PCR product has this length in bp.

# %%
len(prd)

# %% [markdown]
# A figure of the primers annealing on template.

# %%
prd.figure()

# %% [markdown]
# A suggested PCR program.

# %%
prd.program()

# %% [markdown]
# The final vectors are:

# %%
pYPKa_Z_{tp} = (pYPKa_ZraI  + prd).looped().synced("gttctgatcctcgagcatcttaagaattc")
pYPKa_E_{tp} = (pYPKa_EcoRV + prd).looped().synced("gttctgatcctcgagcatcttaagaattc")

# %% [markdown]
# The final vectors with reverse inserts are created below. These vectors theoretically make up
# fifty percent of the clones. The PCR strategy below is used to identify the correct clones.

# %%
pYPKa_Z_{tp}b = (pYPKa_ZraI  + prd.rc()).looped().synced("gttctgatcctcgagcatcttaagaattc")
pYPKa_E_{tp}b = (pYPKa_EcoRV + prd.rc()).looped().synced("gttctgatcctcgagcatcttaagaattc")

# %% [markdown]
# A combination of standard primers and the newly designed primers are
# used for the strategy to identify correct clones.
# Standard primers are listed [here](standard_primers.txt).

# %%
p = {{ x.id: x for x in parse_primers("standard_primers.txt") }}

# %% [markdown]
# ## Diagnostic PCR confirmation
#
# The correct structure of pYPKa_Z_{tp} is confirmed by PCR using standard primers
# 577 and 342 that are vector specific together with the {tp}fw primer specific for the insert
# in a multiplex PCR reaction with
# all three primers present.
#
# Two PCR products are expected if the insert was cloned, the sizes depend
# on the orientation. If the vector is empty or contains another insert, only one
# product is formed.
#
# #### Expected PCR products sizes from pYPKa_Z_{tp}:
#
# pYPKa_Z_{tp} with insert in correct orientation.

# %%
Anneal( (p['577'], p['342'], fp), pYPKa_Z_{tp}).products

# %% [markdown]
# pYPKa_Z_{tp} with insert in reverse orientation.

# %%
Anneal( (p['577'], p['342'], fp), pYPKa_Z_{tp}b).products

# %% [markdown]
# Empty pYPKa clone.

# %%
Anneal( (p['577'], p['342'], fp), pYPKa).products

# %% [markdown]
# #### Expected PCR products sizes pYPKa_E_{tp}:
#
# pYPKa_E_{tp} with insert in correct orientation.

# %%
Anneal( (p['577'], p['342'], fp), pYPKa_E_{tp}).products

# %% [markdown]
# pYPKa_E_{tp} with insert in reverse orientation.

# %%
Anneal( (p['577'], p['342'], fp), pYPKa_E_{tp}b).products


# %% [markdown]
# The cseguid checksum for the resulting plasmids are calculated for future reference.
# The [cseguid checksum](http://pydna.readthedocs.org/en/latest/pydna.html#pydna.utils.cseguid)
# uniquely identifies a circular double stranded sequence.

# %%
print(pYPKa_Z_{tp}.cseguid())
print(pYPKa_E_{tp}.cseguid())

# %% [markdown]
# The sequences are named based on the name of the cloned insert.

# %%
pYPKa_Z_{tp}.locus = "pYPKa_Z_{tp}"[:16]
pYPKa_E_{tp}.locus = "pYPKa_E_{tp}"[:16]

# %% [markdown]
# Sequences are stamped with the cseguid checksum.
# This can be used to verify the integrity of the sequence file.

# %%
pYPKa_Z_{tp}.stamp("cSEGUID")
pYPKa_E_{tp}.stamp("cSEGUID")

# %% [markdown]
# Sequences are written to local files.

# %%
pYPKa_Z_{tp}.write("pYPKa_Z_{tp}.gb")
pYPKa_E_{tp}.write("pYPKa_E_{tp}.gb")
