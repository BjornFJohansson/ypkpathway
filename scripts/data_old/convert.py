#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import notedown, os, nbformat, jupytext
from pathlib import Path

paths = [Path(p) for p in (
            "nb_template_pYPK0_tp_gene_tp.md",
            "nb_template_pYPK0_pw_from_name.md",
            "nb_template_pYPK0_tp_gene_tp_gap_repair.md",
            "nb_template_pYPK0_pw_from_sequence.md",
            "nb_template_pYPKa_ZE_insert.md",
            "nb_template_pYPKa_ZE_genbank_primer_design.md",
            "nb_template_pYPKa_ZE_genbank.md",
            "nb_template_pYPKa_A_insert.md",
            "nb_template_pYPKa_A_genbank_primer_design.md",
            "nb_template_pYPKa_A_genbank.md",
            "nb_template_pYPK0_A_genbank_primer_design.md",
            "nb_template_pYPK0_A_genbank.md")]


obj = notedown.MarkdownReader()

for path in paths:
    newname = path.with_suffix(".py")
    content = path.read_text()   
    jupytext.write(obj.to_notebook(content), str(newname), fmt='py:percent')