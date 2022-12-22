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
# This notebook describes the assembly of {length} transcriptional units
# (single gene expression) vectors into a pathway.
#
# Jupyter notebooks describing the single gene expression vectors are linked
# at the end of this document.
# Specific primers needed are also listed.
#
# ![pathway with N genes](pw.png "pathway with N genes")
# %%
from pydna.parsers import parse_primers
from pydna.readers import read
from pydna.amplify import pcr
from pydna.assembly import Assembly
from IPython.display import display
from IPython.display import Markdown
from pathlib import Path

# %% [markdown]
# The first cassette in the pathway is amplified with standard
# primers 577 and 778, the last with
# 1123 and 578 and all others with 1123 and 778.
# Standard primers are listed [here](standard_primers.fasta).

# %%
p = {{x.name: x for x in parse_primers("standard_primers.fasta")}}

# %% [markdown]
# Restriction enzymes are imported from the Biopython package.

# %%
from Bio.Restriction import {enz}, NotI, PacI

# %% [markdown]
# The backbone vector is linearized by digestion
# with [{enz}](https://www.google.com/search?q={enz}).

# %%
backbone = read("{backbone}.gb")

# %% [markdown]
# The cassette__pcr_products variable holds the list of expression
# cassette PCR products fragments to be assembled.

# %%
cassette_pcr_products = []

# %% [markdown]
# The expression cassettes comes from a series of single gene expression
# vectors held in the template_vectors list.

# %%
cassette_vectors = ("""
{cas_vectors}
""").split()
# %%
cassette_vectors
# %%
template_vectors = [read(v) for v in cassette_vectors]
# %%
for tv in template_vectors:
    display(tv)
# %% [markdown]
# The first cassette in the pathway.
# Suggested PCR conditions can be found at the end of this document.

# %%

fp_first = p['{fp_first}']
fp = p['{fp}']
rp = p['{rp}']
rp_last = p['{rp_last}']

# %%
cassette_pcr_products.append(pcr(fp_first, rp, template_vectors[0]))

# %% [markdown]
# Intermediary cassettes

# %%
cassette_pcr_products.extend(pcr(fp, rp, v)
                             for v in template_vectors[1:-1])

# %% [markdown]
# The last cassette in the pathway.

# %%
cassette_pcr_products.append(pcr(fp, rp_last, template_vectors[-1]))

# %% [markdown]
# The cassettes are given names based on the tu cassette

# %%
for cp, ve in zip(cassette_pcr_products, cassette_vectors):
    cp.name = ve[:-3].split("_", maxsplit=1)[1]
    print(cp.name)

# %% [markdown]
# Cassettes and linear plasmid backbone are joined by homologous recombination

# %%
asm = Assembly([backbone.linearize({enz})] + cassette_pcr_products,
               limit=167-47-10)
asm

# %% [markdown]
# There should normally be two candidates of equal size.
# These sequences should be identical.

# %%
candidates = asm.assemble_circular()
candidates
# %%
candidate, *rest = candidates
# %%
candidate.cseguid() == rest[0].cseguid()
# %% [markdown]
# This assembly figure below shows how the fragments came together.

# %%
candidate.figure()

# %% [markdown]
# The candidate vector is synchronized to the 577 primer. This means that
# the plasmid origin is shifted so that it matches the backbone vector.

# %%
pw = candidate.synced(fp_first)

# %% [markdown]
# The cseguid checksum for the resulting plasmid is calculated for future
# reference.
# The [cseguid checksum](
# http://pydna.readthedocs.org/en/latest/pydna.html#pydna.utils.cseguid)
# uniquely identifies a circular double stranded sequence.

# %%
pw.cseguid()

# %% [markdown]
# The file is given a name based on the sequence of expressed genes.

# %%
pw.locus = "pw"
pw.definition = "{name}"

# %% [markdown]
# Sequence stamped with cseguid checksum.
# This can be used to verify the integrity of the sequence file.

# %%
pw.stamp("cSEGUID")

# %% [markdown]
# Write sequence to a local file.

# %%
pw.write("{name}.gb")

# %% [markdown]
# The pathway can be extended by digestion with either NotI or PacI or both
# provided that the enzymes cut once in the final pathway sequence.

# %%
print(f"NotI cuts {{len(pw.cut(NotI))}} time(s) and PacI cuts "
      f"{{len(pw.cut(PacI))}} time(s) in the final pathway.")

# %% [markdown]
# ### Transcriptional unit (single gene expression) vectors needed.
# %%
for cv in cassette_vectors:
    cassette_vector = Path(cv).with_suffix('.ipynb')
    display(Markdown(f"[{{cassette_vector}}]({{cassette_vector}})"))

# %% [markdown]
# ### Suggested PCR conditions
# %%
for prd in cassette_pcr_products:
    print("\n\n\n\n")
    print("product name:", prd.name)
    print("forward primer", prd.forward_primer.name)
    print("reverse primer", prd.reverse_primer.name)
    print(prd.program())
