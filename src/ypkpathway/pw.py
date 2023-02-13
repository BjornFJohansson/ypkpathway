#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""docstring."""

from pathlib import Path
import re
from pathvalidate import ValidationError, validate_filename

from pydna.dseqrecord import Dseqrecord
from pydna.amplicon import Amplicon
from pydna.primer import Primer

try:
    from importlib.resources import files
except ImportError:  # for Python 3.8
    from importlib_resources import files

from ypkpathway.pth import Pth
from ypkpathway.tu import TranscriptionalUnit

ypkfiles = files("ypkpathway.data")

tu_backbone = "pTA9"


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
