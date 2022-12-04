#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""docstring."""
from importlib.resources import files
import shutil
from pathlib import Path
from typing import Iterable
from datetime import datetime

import jupytext
from nbconvert.preprocessors.execute import ExecutePreprocessor
import nbformat
import pandas

from Bio.Restriction import ZraI, AjiI, EcoRV
from Bio.Restriction import CommOnly as Co

from pydna.amplify import pcr
from pydna.readers import read
from pydna.parsers import parse_primers
from pydna.genbank import genbank


def copy2(src, dst, *, follow_symlinks=True):
    """docstring."""
    try:
        dst = shutil.copy2(src, dst)
    except shutil.SameFileError:
        pass
    return dst

def cloning_primer_design():
    pass




def element_cloning(csv: Path,
                    primerlist: list,
                    cloning_vector: str = "pYPKa",
                    *folders: Iterable[Path],
                    enzymedict={"A": AjiI, "Z": ZraI, "E": EcoRV}):
    """docstring."""

    date = datetime.now().strftime("%d-%b-%Y").upper()

    v = read()  # FIXME: look for sequence file in *folders

    backbonedict = {"Z": v.linearize(ZraI),
                    "A": v.linearize(AjiI),
                    "E": v.linearize(EcoRV)}

    fdict = {"Z": "promoter",
             "A": "gene",
             "E": "terminator"}

    path = Path(csv)
    df = pandas.read_csv(path)
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

        fname = {"Z": f"pYPKa_Z_{row.element}",
                 "A": f"pYPKa_A_{row.element}",
                 "E": f"pYPKa_E_{row.element}"}[row.letter]

        plasmid.locus = fname
        plasmid.definition = fname
        plasmid.accession = "."
        plasmid.annotations["date"] = date
        plasmid.stamp("cSEGUID")
        plasmid.write(f"{fname}.gb")
        paths.append(Path(f"{fname}.gb"))

    return paths


class TranscriptionalUnit:
    """docstring.

    Single-Gene Expression Vector

    class self():pass
    import string;
    set(t[1] for t in string.Formatter().parse(result) if t[1] is not None)
    """

    def __init__(self,
                 name: str = "",
                 *folders: Iterable[Path]):

        self.name = name.removesuffix(".gb")
        self.folders = [Path(f) for f in folders]
        if not self.folders:
            self.folders = [Path.cwd()]
        self.path = None
        self.find_path_to_gb_file()
        self.elements = None
        self.enzymes = None

        self.backbone_name, *elements = self.name.split("_")

        if len(elements) != 3:
            raise ValueError("promoter, gene and terminator needed")

        (self.promoter_name,
         self.gene_name,
         self.terminator_name) = elements

        dp = files("ypkpathway.data")

        self.template_path = dp.joinpath(
                'nb_backbone_promoter_gene_terminator.py')

        self.primers_path = dp.joinpath(
                'standard_primers.fasta')

        self.accessory_file_paths = [dp.joinpath(p) for p in ("tu.png",)]

        assert self.template_path.exists()
        assert self.primers_path.exists()

    def __repr__(self):
        """docstring."""
        return f"tu:{self.name}"

    def find_path_to_gb_file(self):
        """docstring."""
        for folder in self.folders:
            path = (folder/self.name).with_suffix(".gb")
            if path.exists():
                self.path = path
                copy2(path, self.folders[0])
                break
        return self.path

    def find_building_blocks(self):
        """docstring."""
        elements = {self.backbone_name: None,
                    f"pYPKa_Z_{self.promoter_name}": None,
                    f"pYPKa_A_{self.gene_name}": None,
                    f"pYPKa_E_{self.terminator_name}": None}

        for folder in self.folders:
            for element in elements:
                path = (folder/element).with_suffix(".gb")
                if path.exists():
                    elements[element] = path
        return elements

    def copy_files(self):
        """docstring."""
        elements = self.find_building_blocks()
        for name, path in elements.items():
            copy2(path, self.folders[0])
            elements[name] = self.folders[0]/path.name
            assert elements[name].exists()
        return elements

    def check_promoter(self):
        """docstring."""
        p = {x.id: x for x in parse_primers(self.primers_path)}
        promoter_template = read(f"pYPKa_Z_{self.promoter_name}.gb")
        prom = pcr(p['577'], p['567'], promoter_template)
        return prom

    def check_gene(self):
        """docstring."""
        p = {x.id: x for x in parse_primers(self.primers_path)}
        gene_template = read(f"pYPKa_A_{self.gene_name}.gb")
        gene = pcr(p['468'], p['467'], gene_template)
        return gene

    def check_terminator(self):
        """docstring."""
        p = {x.id: x for x in parse_primers(self.primers_path)}
        terminator_template = read(f"pYPKa_E_{self.terminator_name}.gb")
        term = pcr(p['568'], p['578'], terminator_template)
        return term

    def find_enzymes(self):
        """docstring."""
        self.copy_files()        
        bb = read((self.folders[0]/self.backbone_name).with_suffix(".gb"))
        p577, p578, *rest = parse_primers(self.primers_path)
        mcs = pcr(p577, p578, bb)[len(p577):-len(p578)]
        eb = mcs.unique_cutters(Co) & bb.unique_cutters(Co)
        enzymes = [b for a, b in sorted(
            [(abs(len(f1)-len(f2)), e) for ((f1, f2), e) in
             [(mcs.cut(e), e) for e in eb]])]
        return enzymes

    def format_nb(self):
        """docstring."""
        py = self.template_path.read_text().format(
                                            backbone=self.backbone_name,
                                            enz=self.find_enzymes()[0],
                                            promoter=self.promoter_name,
                                            gene=self.gene_name,
                                            terminator=self.terminator_name)

        nb = jupytext.reads(py, fmt='py:percent')

        assert nbformat.validate(nb) is None
        return nb, py

    def execute(self):
        """docstring."""
        self.copy_files()
        copy2(self.primers_path, self.folders[0])
        for p in self.accessory_file_paths:
            copy2(p, self.folders[0])

        nb, py = self.format_nb()

        ep = ExecutePreprocessor()
        ep.timeout = 60  # seconds
        ep.interrupt_on_timeout = True
        resources = {'metadata': {'path': self.folders[0]}}
        nb_executed, resources = ep.preprocess(nb, resources=resources)
        nbformat.write(nb, (self.folders[0]/self.name).with_suffix(".ipynb"))


class PathWay:
    """docstring."""

    def __init__(self,
                 name: str,
                 *folders: Iterable[Path],
                 tu_backbone_name: str = "pTA9"):
        self.name = name
        self.folders = [Path(f) for f in folders]
        self.tu_backbone_name = tu_backbone_name
        self.working_dir = self.folders[0] / Path(name).stem
        dp = files("ypkpathway.data")
        self.template_path = dp.joinpath(
                'nb_backbone_pw_from_name.py')
        self.primers_path = dp.joinpath(
                'standard_primers.fasta')
        self.accessory_file_paths = [dp.joinpath(p) for p in ("pw.png",)]
        self.backbone_name, *elements = name.split("_")

        if len(elements) % 2 != 1:
            raise ValueError("Number of elements must be uneven")

        if len(elements) < 3:  # At least one TU
            raise ValueError("At least a promoter, a gene and a terminator")

        if len(elements) != len(set(e.lower() for e in elements)):
            elms = [e.lower() for e in elements]
            seen = set()
            dupl = []
            for i, e in enumerate(elms):
                if e.lower() in seen:
                    dupl.append(elms.index(e.lower()))
                    dupl.append(i)
                seen.add(e.lower())
            raise ValueError(f"Duplicated elements: "
                             f"{' '.join(elements[i] for i in dupl)}")

        self.transcriptional_units = []

        for i in range(0, len(elements) - 1, 2):
            p, g, t = elements[i:i+3]
            self.transcriptional_units.append(
                TranscriptionalUnit(f"{tu_backbone_name}_{p}_{g}_{t}.gb",
                                    *folders))

        assert (len(self.transcriptional_units) ==
                len([tu.name for tu in self.transcriptional_units]))

    def __repr__(self):
        """docstring."""
        return f"pw:{self.name}"

    def find_enzymes(self):
        """docstring."""
        self.copy_files()
        bb = read((self.folders[0]/self.backbone_name).with_suffix(".gb"))
        p577, p578, *rest = parse_primers(self.primers_path)
        mcs = pcr(p577, p578, bb)[len(p577):-len(p578)]
        eb = mcs.unique_cutters(Co) & bb.unique_cutters(Co)
        enzymes = [b for a, b in sorted(
            [(abs(len(f1)-len(f2)), e) for ((f1, f2), e) in
             [(mcs.cut(e), e) for e in eb]])]
        return enzymes

    def copy_files(self):
        """docstring."""
        elements = self.find_building_blocks()
        for name, path in elements.items():
            copy2(path, self.folders[0])
            elements[name] = self.folders[0]/path.name
            assert elements[name].exists()
        return elements
    
    # def copy_files(self):
    #     """docstring."""
    #     found, missing = self.inspect_project()
    #     for path in found:
    #         copy2(path, self.folders[0])

    def format_nb(self):
        """docstring.

        class self():pass
        import string
        set(t[1] for t in string.Formatter().parse(py) if t[1] is not None)
        """
        cas_vectors = "\n".join(f"{tu.name}.gb" for tu in self.transcriptional_units)
        py = self.template_path.read_text().format(
                                            name=name,
                                            backbone=self.backbone_name,
                                            cas_vectors=cas_vectors,
                                            enz=self.find_enzymes()[0],
                                            length=len(
                                                self.transcriptional_units))

        nb = jupytext.reads(py, fmt='py:percent')

        assert nbformat.validate(nb) is None
        return nb, py

    def execute(self):
        """docstring."""
        for tu in self.transcriptional_units:
            tu.execute()

        copy2(self.primers_path, self.folders[0])
        for p in self.accessory_file_paths:
            copy2(p, self.folders[0])

        nb, py = self.format_nb()

        ep = ExecutePreprocessor()
        ep.timeout = 60  # seconds
        ep.interrupt_on_timeout = True
        resources = {'metadata': {'path': self.folders[0]}}
        nb_executed, resources = ep.preprocess(nb, resources=resources)
        nbformat.write(nb, (self.folders[0]/self.name).with_suffix(".ipynb"))
    

    def find_building_blocks(self):
        """docstring."""
        path = None
        for folder in self.folders:
            path = (folder/self.backbone_name).with_suffix(".gb")
            if path.exists():
                break
        elements = {self.backbone_name: path}
        for tu in self.transcriptional_units:
            tupath = tu.find_path_to_gb_file()
            if tupath:
                elements[tu.name] = tupath
            else:
                tuelements = tu.find_building_blocks()
                for e, p in tuelements.items():
                    if p:
                        elements[e] = p
                    else:
                        elements[e] = None
        return elements


# ipytree https://github.com/QuantStack/ipytree
# anytree https://github.com/c0fec0de/anytree
# treepy https://code.activestate.com/recipes/217212-treepy-graphically-displays-the-directory-structur/


        

if __name__ == "__main__":
    
    # fn = "promoter_list_001.csv"    
    # from pydna.myprimers import PrimerList
    # p = PrimerList()
    # pYPKa_cloning(fn, p)

    folders = ("/home/bjorn/Desktop/ypkpathway/dev/newdev",
               "/home/bjorn/Desktop/ypkpathway/dev")

    name = "pTA1_TDH3_ScCTT1_PGI1"

    name = "pTA1_TDH3_ScATF1_PGI1_ScCTT1_TEF1"

    pw = PathWay(name, *folders)
    
    pw.execute()

    # a,b = pw.inspect_project()
    # NOTE: This is here because sometimes an intermittent issue appears.
    # OPTIMIZE: This could be reworked to not do a O(N2) lookup.
    # TODO: from John: Add a check here to ensure these are always strings.
    # HACK: I am doing something here that is horrible, but it works for now...
    # XXX: Let's do this better next time? It's bad.
    # FIXME: We sometimes get an undefined index in this array.
    # BUG: If the user inputs "Easter" we always output "Egg", even if they wanted a "Bunny".




"""

>568_pCAPsAjiIR (22-mer)
GTGCcatctgtgcagacaaacg

>567_pCAPsAjiIF (23-mer)
GTCggctgcaggtcactagtgag


from anytree import Node, RenderTree
udo = Node("Udo")
marc = Node("Marc", parent=udo)
lian = Node("Lian", parent=marc)
dan = Node("Dan", parent=udo)
jet = Node("Jet", parent=dan)
jan = Node("Jan", parent=dan)
joe = Node("Joe", parent=dan)

for pre, fill, node in RenderTree(udo):
    print("%s%s" % (pre, node.name))
    
    
# https://stackoverflow.com/questions/29850801/subclass-pathlib-path-fails
# (p577,
#  p578,
#  p468,
#  p467,
#  p567,
#  p568,
#  p775,
#  p778,
#  p342) = parse_primers("standard_primers.fasta")





                             >-gene-->
            >-TP-->           \     /           >-TP-->
             \   /             \   /             \   /
     517>     \ /               \ /               \ /
 p577>    1123>|p468>       <p567|p568>       <p467|<494        <p578    <--- recommended
               |                 |                        ep = ExecutePreprocessor()
        ep.timeout = 60  # seconds
        ep.interrupt_on_timeout = Trueb = read(self.backbone_name+".gb")
        resources = {}
        nb_executed, resources = ep.preprocess(nb, resources=resources) |
               |                 |                 |
               |                 |                 |
               |                 |                 |
           775>|                 |                 |<778            
  167>    <511 |<776             |             777>|    <512     <166
               |                 |                 |               <342
               |                 |                 | 
 ✽✽gray✽blue✽N-Z======red========A++++++green+++++E-A••yellow••pink••••••
|            o r                 j                 c c                   |
|            t a                 i                 o c                   |
|            I I                 I                 R I                   |
|                               (*)                V I                   |
|                                                    I                   |
|                                                                        |
 --------------------- pYPKa --------------------------------------------




def inspect_project(self):
    found = []
    missing = []
    for tu in self.transcriptional_units:
        tupath  = tu.find_path_to_gb_file()
        if tupath:
            found.append(tupath)
        else:
            elements = tu.find_building_blocks()
            for e, p in elements.items():
                if p:
                    found.append(p)
                else:
                    missing.append(p)
    return found, missing

"""











