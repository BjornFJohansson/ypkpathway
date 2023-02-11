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
# # Construction of promoter vector pYPKa_Z_{insert} and terminator vector pYPKa_E_{insert}
#
# This notebook describe the construction of _E. coli_ vectors [pYPKa_Z_{insert}](pYPKa_Z_{insert}.gb) and [pYPKa_E_{insert}](pYPKa_E_{insert}.gb)
# with the same insert for which PCR primers are also designed.
#
# The insert defined below is cloned in pYPKa using blunt restriction
# enzymes [ZraI](http://rebase.neb.com/rebase/enz/ZraI.html) or [EcoRV](http://rebase.neb.com/rebase/enz/EcoRV.html) resulting in
# two different plasmids. The insert cloned in ZraI will be used as a promoter, while the insert in EcoRV will be used as a terminator.
#
# ![pYPKa_Z and pYPKa_E plasmids](pYPK_ZE.png "pYPKa_Z and pYPKa_E plasmids")
#
# A part of the [pydna](https://pypi.python.org/pypi/pydna/) package is imported in the code cell below.

# %%
from pydna.readers import read
from pydna.genbank import Genbank
from pydna.parsers import parse_primers
from pydna.design  import primer_design
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
# Primers are needed to PCR amplify the insert. The forward primer adds two adenines in front of the start codon
# which is a feature commonly found in highly expressed _S. cerevisiae_ genes.

# %%
fp_tail = "ttaaat"
rp_tail = "taattaa"

# %% [markdown]
# Primers with the tails above are designed in the code cell below.

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
# The PCR product has this length in bp.

# %%
len(ins)

# %% [markdown]
# The primers anneal on the template like this.

# %%
ins.figure()

# %% [markdown]
# A suggested PCR program.

# %%
ins.program()

# %% [markdown]
# The final vectors are:

# %%
pYPKa_Z_{insert} = (pYPKa_ZraI  + ins).looped().synced("gttctgatcctcgagcatcttaagaattc")
pYPKa_E_{insert} = (pYPKa_EcoRV + ins).looped().synced("gttctgatcctcgagcatcttaagaattc")

# %% [markdown]
# The final vectors with reverse inserts are created below. These vectors theoretically make up
# fifty percent of the clones. The PCR strategy below is used to identify the correct clones.

# %%
pYPKa_Z_{insert}b = (pYPKa_ZraI  + ins.rc()).looped().synced("gttctgatcctcgagcatcttaagaattc")
pYPKa_E_{insert}b = (pYPKa_EcoRV + ins.rc()).looped().synced("gttctgatcctcgagcatcttaagaattc")

# %% [markdown]
# A combination of standard primers and the newly designed primers are
# used for the strategy to identify correct clones.
# Standard primers are listed [here](standard_primers.txt).

# %%
p = {{ x.id: x for x in parse_primers("standard_primers.txt") }}

# %% [markdown]
# ## Diagnostic PCR confirmation
#
# The correct structure of pYPKa_Z_{insert} is confirmed by PCR using standard primers
# 577 and 342 that are vector specific together with the {insert}fw primer specific for the insert
# in a multiplex PCR reaction with
# all three primers present.
#
# Two PCR products are expected if the insert was cloned, the sizes depend
# on the orientation. If the vector is empty or contains another insert, only one
# product is formed.
#
# #### Expected PCR products sizes from pYPKa_Z_{insert}:
#
# pYPKa_Z_{insert} with insert in correct orientation.

# %%
Anneal( (p['577'], p['342'], fp), pYPKa_Z_{insert}).products

# %% [markdown]
# pYPKa_Z_{insert} with insert in reverse orientation.

# %%
Anneal( (p['577'], p['342'], fp), pYPKa_Z_{insert}b).products

# %% [markdown]
# Empty pYPKa clone.

# %%
Anneal( (p['577'], p['342'], fp), pYPKa).products

# %% [markdown]
# #### Expected PCR products sizes pYPKa_E_{insert}:
#
# pYPKa_E_{insert} with insert in correct orientation.

# %%
Anneal( (p['577'], p['342'], fp), pYPKa_E_{insert}).products

# %% [markdown]
# pYPKa_E_{insert} with insert in reverse orientation.

# %%
Anneal( (p['577'], p['342'], fp), pYPKa_E_{insert}b).products

# %% [markdown]
# The cseguid checksum for the resulting plasmids are calculated for future reference.
# The [cseguid checksum](http://pydna.readthedocs.org/en/latest/pydna.html#pydna.utils.cseguid)
# uniquely identifies a circular double stranded sequence.

# %%
print(pYPKa_Z_{insert}.cseguid())
print(pYPKa_E_{insert}.cseguid())

# %% [markdown]
# The sequences are named based on the name of the cloned insert.

# %%
pYPKa_Z_{insert}.locus = "pYPKa_Z_{insert}"[:16]
pYPKa_E_{insert}.locus = "pYPKa_E_{insert}"[:16]

# %% [markdown]
# Sequences are stamped with the cseguid checksum.
# This can be used to verify the integrity of the sequence file.

# %%
pYPKa_Z_{insert}.stamp("cSEGUID")
pYPKa_E_{insert}.stamp("cSEGUID")

# %% [markdown]
# Sequences are written to local files.

# %%
pYPKa_Z_{insert}.write("pYPKa_Z_{insert}.gb")
pYPKa_E_{insert}.write("pYPKa_E_{insert}.gb")
