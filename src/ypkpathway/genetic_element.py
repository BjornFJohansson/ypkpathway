#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""docstring."""

from pathlib import Path
from typing import Iterable
from datetime import datetime

import pandas

from Bio.Restriction import ZraI, AjiI, EcoRV

from Bio.SeqRecord import SeqRecord

from pydna.amplify import pcr
from pydna.readers import read
from pydna.genbank import genbank


def _find_path_to_file(folders: Iterable[Path], fn: str):
    """docstring."""
    foundpath = None
    for folder in folders:
        path = folder/fn
        if path.exists():
            foundpath = path
            break
    return foundpath


def element_cloning(csv: Path,
                    datafolders: Iterable[Path] = (Path("."),),
                    primerlist: Iterable[SeqRecord] = tuple(),
                    cloning_vector: str = "pYPKa",
                    enzymedict={"A": AjiI, "Z": ZraI, "E": EcoRV}):
    """docstring."""
    date = datetime.now().strftime("%d-%b-%Y").upper()

    datafolders = [Path(f) for f in datafolders]

    cloning_vector = str(Path(cloning_vector).stem)

    vpth = _find_path_to_file(datafolders, cloning_vector)

    v = read(vpth)

    backbonedict = {"Z": v.linearize(ZraI),
                    "A": v.linearize(AjiI),
                    "E": v.linearize(EcoRV)}

    fdict = {"Z": "promoter",
             "A": "gene",
             "E": "terminator"}

    csvpath = Path(csv)
    df = pandas.read_csv(csvpath)
    paths = []
    p = primerlist

    for i, row in df.iterrows():
        insert = pcr(p[row.fp], p[row.rp], genbank(row.genbank))
        label = f"{row.element}_{fdict[row.letter]}"

        insert.add_feature(type="gene",
                           name=label,
                           ApEinfo_label=label,
                           ApEinfo_fwdcolor="#FC6600",
                           ApEinfo_revcolor="#7C4700",
                           fp=p[row.fp].format("tab").strip(),
                           rp=p[row.rp].format("tab").strip(),
                           reference=row.genbank)
        if row.letter != "A":
            f = insert.features[-1]
            f.type = "regulatory"
            qf = {"regulatory_class": fdict[row.letter]}
            qf.update(f.qualifiers)
            f.qualifiers = qf

        plasmid = (backbonedict[row.letter] + insert).looped().synced(p[577])

        fname = f"pYPKa_{row.letter.upper()}_{row.element}",

        plasmid.locus = fname
        plasmid.definition = fname
        plasmid.accession = "."
        plasmid.annotations["date"] = date
        plasmid.stamp("cSEGUID")
        plasmid.write(f"{fname}.gb")
        paths.append(Path(f"{fname}.gb"))

    return paths


if __name__ == "__main__":

    fn = "promoter_list_001.csv"
    from pydna.myprimers import PrimerList
    p = PrimerList()
    element_cloning(fn, p)
