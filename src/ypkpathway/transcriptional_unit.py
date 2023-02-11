#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""docstring."""

from pathlib import Path
from pprint import pprint
import re
from pathvalidate import ValidationError, validate_filename
# from copy import deepcopy

import jupytext
from nbconvert.preprocessors.execute import ExecutePreprocessor
import nbformat
from Bio.Restriction import CommOnly as Co
from pydna.amplify import pcr
from pydna.readers import read
from pydna.parsers import parse_primers
from pydna.assembly import Assembly

from ypkpathway.pth import Pth

from pydna.dseqrecord import Dseqrecord
from pydna.amplicon import Amplicon
from pydna.primer import Primer

try:
    from importlib.resources import files
except ImportError:  # for Python 3.8
    from importlib_resources import files

ypkfiles = files("ypkpathway.data")

tu_backbone = "pTA9"


class TranscriptionalUnit:
    """docstring."""

    def __init__(self,
                 name: str,
                 workdir: Path = Path("."),
                 datadirs: list[Path] = None,
                 structure: dict = None,
                 *args,
                 **kwargs):

        tu = {"path": Pth("name"),
              "backbone": Dseqrecord("", id="na", name="backbone"),
              "parts": [Amplicon("",
                                 id="pYPKa_Z_",
                                 name="promoter",
                                 forward_primer=Primer("",
                                                       name="577_crp585-557"),
                                 reverse_primer=Primer("",
                                                       name="567_pCAPsAjiIF")),
                        Amplicon("",
                                 id="pYPKa_A_",
                                 name="gene",
                                 forward_primer=Primer("",
                                                       name="468_pCAPs_release_fw"),
                                 reverse_primer=Primer("",
                                                       name="467_pCAPs_release_re")),
                        Amplicon("",
                                 id="pYPKa_E_",
                                 name="terminator",
                                 forward_primer=Primer("",
                                                       name="568_pCAPsAjiIR"),
                                 reverse_primer=Primer("",
                                                       name="578_crp42-70")), ],
              "primer_file_path": Pth("standard_primers.fasta", "."),
              "template_file_path": ypkfiles/"nb_tu.py",
              "accessory_file_paths": [Pth(p, (ypkfiles,)) for p in ("tu.png",
                                                                     "logo.png")]}

        structure = structure or tu
        backbone, elements = self._generate_elements(name)
        self.__dict__.update(structure)
        datadirs = datadirs or [Path("."), ]
        datadirs.insert(0, workdir)
        datadirs = [Path(d) for d in datadirs]
        workdir = Pth(workdir, datadirs)/name
        self.path = (workdir/name).with_suffix(".gb")
        self.backbone.id = backbone
        self.path.parent.mkdir(parents=True, exist_ok=True)
        if list(self.path.parent.glob("*")):
            raise Exception(f"{workdir} exists and is not empty.")
        assert all(d.is_dir() for d in self.path.datadirs)
        self.primer_file_path = workdir/self.primer_file_path
        # self.template_file_path = str(workdir)/self.template_file_path
        for i, p in enumerate(self.accessory_file_paths):
            self.accessory_file_paths[i] = str(workdir)/p
        # self.parts = deepcopy(self.parts)
        # self.accessory_file_paths = deepcopy(self.accessory_file_paths)
        for elementname, element in zip(elements, self.parts):
            # self.parts[i] = deepcopy(element)
            # self.parts[i].id += elementname
            element.id += elementname

    def _generate_elements(self, name):
        try:
            # The name has to be a valid filename
            validate_filename(name, check_reserved=True)
        except ValidationError as e:
            raise e
        if "." in name:
            raise Exception("dot (.) not permitted in name")
        backbone, *elements = name.split("_")
        if len(elements) != 3:
            raise ValueError("Name has to indicate a backbone, a promoter "
                             "a gene and a terminator; four elements "
                             "separated by underscore (_). "
                             f"Not: {name}")
        return backbone, elements

    def collect(self):
        """docstring."""
        primer_path = self.primer_file_path
        primer_path.copy_to_dir()
        pd = {x.name: x for x in parse_primers(primer_path)}

        for i, element in enumerate(self.parts):
            name = element.id
            p = Pth((self.path.parent/name).with_suffix(".gb"),
                    self.path.datadirs)
            p.copy_to_dir()
            try:
                element.template = read(p)
            except ValueError:
                print(f"no {p} found or no sequence")
            else:
                element.forward_primer = pd.get(element.forward_primer.name)
                element.reverse_primer = pd.get(element.reverse_primer.name)
                self.parts[i] = element

        p = (self.path.parent/self.backbone.id).with_suffix(".gb")
        p.copy_to_dir()
        bbs = read(p)
        f = self.parts[0].forward_primer
        r = self.parts[-1].reverse_primer
        mcs = pcr(f, r, bbs)[len(f):-len(r)]
        eb = mcs.unique_cutters(Co) & bbs.unique_cutters(Co)
        enzymes = [b for a, b in sorted(
              [(abs(len(f1)-len(f2)), e) for ((f1, f2), e) in
               [(mcs.cut(e), e) for e in eb]])]
        bbs.annotations["comment"] = f"{bbs.annotations.get('comment', '')}\nLinearized with {enzymes[0]}"
        bbs.id = self.backbone.id
        self.backbone = bbs

        for obj in self.__dict__.values():
            if hasattr(obj, "copy_to_dir"):
                obj.copy_to_dir()
            elif isinstance(obj, list):
                for elm in obj:
                    if hasattr(elm, "copy_to_dir"):
                        elm.copy_to_dir()

    def _format_code(self):
        """docstring."""
        py = self.template_file_path.read_text()
        *rest, last = self.backbone.annotations["comment"].splitlines()
        *rest, enz = last.split()

        pycode = py.format(name=repr(self),
                           backbone=self.backbone.id,
                           enz=enz,
                           promoter=self.parts[0].id,
                           gene=self.parts[1].id,
                           terminator=self.parts[2].id,
                           fp_prom=self.parts[0].forward_primer.name,
                           rp_prom=self.parts[0].reverse_primer.name,
                           fp_gene=self.parts[1].forward_primer.name,
                           rp_gene=self.parts[1].reverse_primer.name,
                           fp_term=self.parts[2].forward_primer.name,
                           rp_term=self.parts[2].reverse_primer.name)
        
        return pycode
    
    def generate_notebook(self):
        py = self._format_code()
        nb = jupytext.reads(py, fmt='py:percent')
        assert nbformat.validate(nb) is None
        ep = ExecutePreprocessor()
        ep.timeout = 60  # seconds
        ep.interrupt_on_timeout = True
        output = ""
        resources = {'metadata': {'path': self.path.parent, 
                                  "stdout": output}}
        nb_executed, resources = ep.preprocess(nb, 
                                               resources=resources)
        nbformat.write(nb, self.path.with_suffix(".ipynb"))
        return resources

    def __repr__(self):
        """docstring."""
        return str(self.path.stem)

    @property
    def status(self):
        """docstring."""
        pprint(self.__dict__, sort_dicts=False)

    def assemble(self):
        """docstring."""
        prom = pcr(self.parts[0])
        gene = pcr(self.parts[1])
        term = pcr(self.parts[2])
        *rest, last = self.backbone.annotations["comment"].splitlines()
        *rest, enz = last.split()
        from Bio.Restriction import RestrictionBatch
        linear_backbone = self.backbone.linearize(RestrictionBatch([enz]))
        asm = Assembly((linear_backbone, prom, gene, term), limit=31)
        candidates = asm.assemble_circular()
        candidate, *rest = candidates
        return candidate, rest


class PathWay(TranscriptionalUnit):
    """docstring."""

    def __init__(self,
                 name: str,
                 workdir: Path = Path("."),
                 datadirs: list[Path] = None,
                 structure: dict = None,
                 *args,
                 **kwargs):
        
        pw = {"path": Pth("name"),
              "backbone": Dseqrecord("", id="na", name="backbone"),
              "parts": [Amplicon("",
                                 id=f"{tu_backbone}_",
                                 name="first_cassette",
                                 forward_primer=Primer("",
                                                       name="577_crp585-557"),
                                 reverse_primer=Primer("",
                                                       name="778_tp_Eco32I_rev")),
                        Amplicon("",
                           id=f"{tu_backbone}_",
                           name="middle_cassette",
                           forward_primer=Primer("",
                                                 name="1123_New775"),
                           reverse_primer=Primer("",
                                                 name="778_tp_Eco32I_rev")),
                        Amplicon("",
                                 id=f"{tu_backbone}_",
                                 name="last_cassette",
                                 forward_primer=Primer("",
                                                       name="1123_New775"),
                                 reverse_primer=Primer("",
                                                       name="578_crp42-70")), ],
              "primer_file_path": Pth("standard_primers.fasta", "."),
              "template_file_path": ypkfiles/"nb_pw.py",
              "accessory_file_paths": [Pth(p, (ypkfiles,)) for p in ("pw.png",
                                                                     "logo.png")]}
        structure = structure or pw
        super().__init__(name,
                         workdir,
                         datadirs,
                         structure,
                         *args, 
                         **kwargs)
        
    def collect(self):
        """docstring."""
        tunits = []
        for part in self.parts:
            p = (self.path.parent/part.id).with_suffix(".gb")
            if p.find_in_dirs():
                p.copy_to_dir()
            else:
                tunits.append(TranscriptionalUnit(part.id,
                                                  self.path.parent,
                                                  self.path.datadirs))
        for tu in tunits:
            tu.collect()
            tu.generate_notebook()

        super().collect()
            


    def _generate_elements(self, name):
        try:
            # The name has to be a valid filename
            validate_filename(name, check_reserved=True)
        except ValidationError as e:
            raise e
        if "." in name:
            raise Exception("dot (.) not permitted in name")

        backbone, *elements = name.split("_")

        if len(elements) < 5 or len(elements) % 2 != 1:
            raise ValueError("At least five uneven number of elements "
                             "(two TU casettes) separated by underscore (_). "
                             f"Not {name}")

        if len(elements) != len(set(e.lower() for e in elements)):
            elms = [e.lower() for e in elements]
            seen = set()
            dupl = []
            for i, e in enumerate(elms):
                if e.lower() in seen:
                    dupl.append(elms.index(e.lower()))
                    dupl.append(i)
                seen.add(e.lower())
            pattern = name
            for elm in set(elements[i] for i in dupl):
                pattern = pattern.replace(elm, len(elm)*"*")
            pattern = re.sub("[^*]", " ", pattern)
            raise ValueError(f"Duplicated elements: "
                             f"{' '.join(set(elements[i] for i in dupl))}\n"
                             f"{name}\n"
                             f"{pattern}\n")
        newelements = []
        for i in range(0, len(elements) - 1, 2):
            n = "{}_{}_{}".format(*elements[i:i+3])
            newelements.append(n)
        return backbone, newelements

    def _format_code(self):
        """docstring."""
        py = self.template_file_path.read_text()
        *rest, last = self.backbone.annotations["comment"].splitlines()
        *rest, enz = last.split()
        pycode = py.format(name=repr(self),
                           backbone=self.backbone.id,
                           cas_vectors="\n".join(f"{p.id}.gb" for p in self.parts),
                           enz=enz,
                           length=len(self.parts),
                           fp_first=self.parts[0].forward_primer.name,
                           fp=self.parts[-1].forward_primer.name,
                           rp=self.parts[0].reverse_primer.name,
                           rp_last=self.parts[-1].reverse_primer.name)        
        return pycode

if __name__ == "__main__":
    
    import shutil

    shutil.rmtree("/home/bjorn/python_packages/ypkpathway/src/ypkpathway/pYPK0_PDC1_KlLAC4_PGI1_KlLAC12_TPI1_ScELO1_TDH3", ignore_errors=True)

    workdir = "/home/bjorn/Desktop/python_packages/ypkpathway/src/ypkpathway"

    datadirs = [
        "/home/bjorn/Desktop/mec@githb/YeastPathwayKit/sequences",
        "/home/bjorn/Desktop/mec@githb/YeastPathwayKitPrivate/sequences",]

    name = "pYPK0_PDC1_KlLAC4_PGI1"
    name = "pYPK0_PDC1_KlLAC4_PGI1_KlLAC12_TPI1"
    name = "pYPK0_PDC1_KlLAC4_PGI1_KlLAC12_TPI1_ScELO1_TDH3"
    
    self = PathWay(name, workdir, datadirs)
    self.collect()
    self.collect()
    self.generate_notebook()
    
    
    # t = TranscriptionalUnit("pTA9_PDC1_KlLAC4_PGI1", workdir, datadirs)
    
    # [x.id for x in t.parts]
    
    
    
    # self.collect()
    # self.assemble()
    # self.generate_notebook()
