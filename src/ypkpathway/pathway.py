#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""docstring."""
import shutil
import re
from pathlib import Path
from typing import Iterable
from datetime import datetime as dt
from pprint import pprint

import jupytext
from nbconvert.preprocessors.execute import ExecutePreprocessor
import nbformat

from pathvalidate import ValidationError, validate_filename

from Bio.Restriction import CommOnly as Co

from pydna.amplify import pcr
from pydna.readers import read
from pydna.parsers import parse_primers

try:
    from importlib.resources import files
except ImportError:  # for Python 3.8
    from importlib_resources import files

from ypkpathway.transcriptional_unit import TranscriptionalUnit
from ypkpathway.utils import _copy2, _find_path_to_file


class PathWay:
    """docstring."""

    def __init__(self,
                 name: str,
                 workdir: Path = Path("."),
                 datafolders: Iterable[Path] = (Path("."),),
                 tu_backbone: str = "pTA9",
                 primernames={"fp_first": "577_crp585-557",
                              "fp": "1123_New775",
                              "rp": "778_tp_Eco32I_rev",
                              "rp_last": "578_crp42-70", },
                 primersfn="standard_primers.fasta"):

        try:
            validate_filename(name)
        except ValidationError as e:
            raise e

        if not name:
            raise ValueError("name is mandatory.")

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

        self.paths = {"fn": Path(name).with_suffix(".gb")}
        self.paths["backbone"] = Path(backbone).with_suffix(".gb")

        self.transcriptional_units = {}

        workdir = Path(workdir)/name

        for i in range(0, len(elements) - 1, 2):
            n = "{}_{}_{}_{}".format(tu_backbone, *elements[i:i+3])
            self.transcriptional_units[n] = TranscriptionalUnit(
                                                n,
                                                workdir,
                                                datafolders)

        self.paths.update({n: tu.paths["fn"]
                           for n, tu in self.transcriptional_units.items()})

        self.files = {key: None for key in list(self.paths.keys())}

        self.products = dict(list(self.files.items())[2:])

        self.primers = {}
        self.primers["names"] = primernames
        self.primers["fn"] = Path(primersfn)

        self.workdir = workdir
        self.datafolders = [self.workdir] + [Path(f) for f in datafolders]

        self.enzyme = None
        self.workdir.mkdir(parents=True, exist_ok=True)

        dp = files("ypkpathway.data")
        self.template_path = dp.joinpath('nb_backbone_pw_from_name.py')
        self.accessory_file_paths = [dp.joinpath(p) for p in ("pw.png",
                                                              "logo.png")]

    def __repr__(self):
        """docstring."""
        return str(self.paths["fn"].stem)

    @property
    def status(self):
        """docstring."""
        pprint(self.__dict__, sort_dicts=False)

    def find_sequence_files(self):
        """docstring."""
        notfound = {}
        for key, fn in self.paths.items():
            path = _find_path_to_file(self.datafolders, fn)
            if path:
                self.paths[key] = path
            else:
                notfound[key] = fn
        return notfound

    def copy_sequence_files(self):
        """docstring."""
        not_copied = {}
        for key, path in self.paths.items():
            if path.is_file():
                p = _copy2(path, self.workdir)
                # breakpoint()
                # assert p.is_file()
            else:
                not_copied[key] = path
        return not_copied

    def read_sequence_files(self):
        """docstring."""
        not_read = {}
        for element, path in self.paths.items():
            try:
                self.files[element] = read(path)
            except ValueError as err:
                if err.args[0].startswith("No sequences found in data:"):
                    not_read[element] = path
                    pass
            except Exception as err:
                print(f"Unexpected {err=}, {type(err)=}")
                raise
        return not_read

    def find_primers(self):
        """docstring."""
        path = _find_path_to_file(self.datafolders, self.primers["fn"])
        if path:
            self.primers["fn"] = path
            path = None
        return path

    def copy_primers(self):
        """docstring."""
        path = _copy2(self.primers["fn"], self.workdir)
        if path.is_file():
            self.primers["fn"] = path
            msg = None
        else:
            msg = path
        return msg

    def read_primers(self):
        """docstring."""
        self.primers["fn"] = _find_path_to_file(self.datafolders,
                                                self.primers["fn"])
        try:
            primerdict = {x.name: x
                          for x in
                          parse_primers(self.primers["fn"])}
        except Exception as err:
            print(f"Unexpected {err=}, {type(err)=}")
            raise

        for name in self.primers["names"].values():
            self.primers[name] = primerdict.get(name, None)

    def find_enzymes(self):
        """docstring."""
        bb = self.files["backbone"]
        f = self.primers[self.primers["names"]["fp_first"]]
        r = self.primers[self.primers["names"]["rp_last"]]
        mcs = pcr(f, r, bb)[len(f):-len(r)]
        eb = mcs.unique_cutters(Co) & bb.unique_cutters(Co)
        enzymes = [b for a, b in sorted(
            [(abs(len(f1)-len(f2)), e) for ((f1, f2), e) in
             [(mcs.cut(e), e) for e in eb]])]
        self.enzyme = enzymes[0]
        return enzymes

    def execute_transcriptional_unit_notebooks(self):
        """docstring."""
        from tqdm import tqdm
        for tuname, tu in tqdm(self.transcriptional_units.items()):
            tu.find_primers()
            tu.copy_primers()
            tu.read_primers()
            tu.find_sequence_files()
            tu.copy_sequence_files()
            tu.read_sequence_files()
            tu.find_enzymes()
            tu.execute_nb()
            self.paths[tuname] = tu.paths["fn"]

    def execute_transcriptional_unit_python(self):
        """docstring."""
        pass

    def pcr(self):
        """docstring."""
        f = self.primers[self.primers["names"]["fp_first"]]
        r = self.primers[self.primers["names"]["rp"]]
        first_cassette = list(self.files.values())[3]
        self.products[list(self.files.keys())[3]] = pcr(f, r, first_cassette)

        f = self.primers[self.primers["names"]["fp"]]
        r = self.primers[self.primers["names"]["rp_last"]]
        last_cassette = list(self.files.values())[-1]
        self.products[list(self.files.keys())[-1]] = pcr(f, r, last_cassette)

        f = self.primers[self.primers["names"]["fp"]]
        r = self.primers[self.primers["names"]["rp"]]

        for key, template in list(self.files.items())[3:-1]:
            cassette = list(self.files.values())[-1]
            self.products[key] = pcr(f, r, cassette)

    def format_code(self):
        """docstring."""
        cas_vectors = "\n".join(tu.paths["fn"].name for tu in
                                self.transcriptional_units.values())
        py = self.template_path.read_text().format(
                                            name=repr(self),
                                            backbone=self.paths["backbone"].stem,
                                            cas_vectors=cas_vectors,
                                            enz=self.enzyme,
                                            length=len(
                                                self.transcriptional_units),
                                            fp_first=self.primers["names"]['fp_first'],
                                            fp=self.primers["names"]['fp'],
                                            rp=self.primers["names"]['rp'],
                                            rp_last=self.primers["names"]['rp_last'])

        nb = jupytext.reads(py, fmt='py:percent')

        assert nbformat.validate(nb) is None
        return nb, py

    def execute_py(self):
        """docstring."""
        for p in self.accessory_file_paths:
            _copy2(p, self.workdir)
        self.copy_primers()
        self.copy_sequence_files()
        globs = {}
        nb, py = self.format_code()
        cwd = Path.cwd()
        os.chdir(self.workdir)
        f = StringIO()
        with redirect_stdout(f):
            exec(py, globs)
        s = f.getvalue()        
        os.chdir(cwd)
        name = Path(self.workdir/self.paths["fn"]).with_suffix(".py").write_text(py)
        return s

    def execute_nb(self):
        """docstring."""
        for p in self.accessory_file_paths:
            _copy2(p, self.workdir)
        self.copy_primers()
        self.copy_sequence_files()
        nb, py = self.format_code()
        ep = ExecutePreprocessor()
        ep.timeout = 60  # seconds
        ep.interrupt_on_timeout = True
        output = ""
        resources = {'metadata': {'path': self.workdir, 
                                  "stdout": output}}
        nb_executed, resources = ep.preprocess(nb, 
                                               resources=resources)
        
        nbformat.write(nb, (self.workdir/self.paths["fn"]).with_suffix(".ipynb"))
        
        return resources

    def assemble(self):
        prom = self.pcr_products["promoter"]
        gene = self.pcr_products["gene"]
        term = self.pcr_products["terminator"]

        prom.name = self.paths["promoter"].stem[-4:]
        gene.name = self.paths["gene"].stem[-4:]
        term.name = self.paths["terminator"].stem[-4:]

        linear_backbone = self.files["backbone"].linearize(self.enzyme)

        asm = Assembly((linear_backbone, prom, gene, term), limit=31)

        candidates = asm.assemble_circular()

        candidate, *rest = candidates

        candidate.cseguid() == rest[0].cseguid()

        fp_prom = self.primers[self.primernames["fp_prom"]]

        return candidate.synced(fp_prom), candidate.figure()

# ipytree https://github.com/QuantStack/ipytree
# anytree https://github.com/c0fec0de/anytree
# treepy https://code.activestate.com/recipes/217212-treepy-graphically-displays-the-directory-structur/

if __name__ == "__main__":

    workdir = "/home/bjorn/Desktop/python_packages/ypkpathway/tests/test/new"

    datafolders = (
        "/home/bjorn/Desktop/python_packages/ypkpathway/tests/test",
        "/home/bjorn/Desktop/YeastPathwayKit/sequences",)

    name = "pTA1_TDH3_ScATF1_PGI1_ScCTT1_TEF1"
    
    # name = "pTA1_TDH3_ScATF1_PGI1_ScCTT1_TEF1_TDH3_pgi1"

    self = PathWay(name, 
                   workdir, 
                   datafolders)

    self.find_sequence_files()    
    self.copy_sequence_files()    
    self.read_sequence_files()    
    self.find_primers()
    self.copy_primers()
    self.read_primers()    
    self.find_enzymes()    
    # self.execute_transcriptional_unit_notebooks()    
    self.find_sequence_files()    
    self.read_sequence_files()    
    self.pcr()
    nb, py = self.format_code()
    print(py)
    
    self.execute_nb()
