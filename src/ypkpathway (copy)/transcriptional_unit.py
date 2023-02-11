#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""docstring."""

from pathlib import Path
from typing import Iterable
from pprint import pprint
import os
from io import StringIO
from contextlib import redirect_stdout

import jupytext
from nbconvert.preprocessors.execute import ExecutePreprocessor
import nbformat

from Bio.Restriction import CommOnly as Co

from pydna.amplify import pcr
from pydna.readers import read
from pydna.parsers import parse_primers
from pydna.assembly import Assembly

from pathvalidate import ValidationError, validate_filename

try:
    from importlib.resources import files
except ImportError:  # for Python 3.8
    from importlib_resources import files

from ypkpathway.utils import _copy2, _find_path_to_file


class TranscriptionalUnit:
    """docstring.

    Single-Gene Expression Vector

    """

    def __init__(self,
                 name: str = "",
                 workdir: Path = Path("."),
                 datafolders: Iterable[Path] = (Path("."),),
                 primernames={"fp_prom": "577_crp585-557",
                              "rp_prom": "567_pCAPsAjiIF",
                              "fp_gene": "468_pCAPs_release_fw",
                              "rp_gene": "467_pCAPs_release_re",
                              "fp_term": "568_pCAPsAjiIR",
                              "rp_term": "578_crp42-70", },
                 primersfn="standard_primers.fasta"):

        try:
            validate_filename(name)
        except ValidationError as e:
            raise e

        backbone, *elements = name.split("_")

        if len(elements) != 3:
            raise ValueError("Name has to indicate a backbone, a promoter "
                             "a gene and a terminator; four elements "
                             "separated by underscore (_). "
                             f"Not: {name}")

        self.paths = {"fn": Path(name).with_suffix(".gb")}
        self.paths["backbone"] = Path(backbone).with_suffix(".gb")
        self.paths["promoter"] = Path(f"pYPKa_Z_{elements[0]}").with_suffix(".gb")
        self.paths["gene"] = Path(f"pYPKa_A_{elements[1]}").with_suffix(".gb")
        self.paths["terminator"] = Path(f"pYPKa_E_{elements[2]}").with_suffix(".gb")

        self.files = {key: None for key in list(self.paths.keys())}

        self.pcr_products = dict(list(self.files.items())[2:])

        self.primers = {}
        self.primers["names"] = primernames
        self.primers["fn"] = Path(primersfn)

        self.workdir = Path(workdir)
        self.datafolders = [Path(f) for f in datafolders]

        self.enzyme = None
        self.workdir.mkdir(parents=True, exist_ok=True)

        dp = files("ypkpathway.data")
        self.template_path = dp.joinpath(
                'nb_backbone_promoter_gene_terminator.py')
        assert self.template_path.exists()
        self.accessory_file_paths = [dp.joinpath(p) for p in ("tu.png",
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
                assert p.is_file()
            else:
                not_copied[key] = path
        return not_copied

    def read_sequence_files(self):
        """docstring."""
        for element, path in self.paths.items():
            try:
                self.files[element] = read(path)
            except ValueError as err:
                if err.args[0].startswith("No sequences found in data:"):
                    pass
            except Exception as err:
                print(f"Unexpected {err=}, {type(err)=}")
                raise

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
        f = self.primers[self.primers["names"]["fp_prom"]]
        r = self.primers[self.primers["names"]["rp_term"]]
        mcs = pcr(f, r, bb)[len(f):-len(r)]
        eb = mcs.unique_cutters(Co) & bb.unique_cutters(Co)
        enzymes = [b for a, b in sorted(
            [(abs(len(f1)-len(f2)), e) for ((f1, f2), e) in
             [(mcs.cut(e), e) for e in eb]])]
        self.enzyme = enzymes[0]
        return enzymes

    def pcr(self):
        """docstring."""
        
        fp_prom = self.primers[self.primers["names"]["fp_prom"]]
        rp_prom = self.primers[self.primers["names"]["rp_prom"]]
        
        self.pcr_products["promoter"] = pcr(fp_prom,
                                            rp_prom,
                                            self.files['promoter'])
        
        fp_gene = self.primers[self.primers["names"]["fp_gene"]]
        rp_gene = self.primers[self.primers["names"]["rp_gene"]]
        
        self.pcr_products["gene"] = pcr(fp_gene,
                                        rp_gene,
                                        self.files['gene'])
        
        fp_term = self.primers[self.primers["names"]["fp_term"]]
        rp_term = self.primers[self.primers["names"]["rp_term"]]    

        self.pcr_products["terminator"] = pcr(fp_term,
                                              rp_term,
                                              self.files['terminator'])

    def format_code(self):
        """docstring."""
        py = self.template_path.read_text().format(
                                            name=repr(self),
                                            backbone=self.paths["backbone"].stem,
                                            enz=self.enzyme,
                                            promoter=self.paths["promoter"].stem,
                                            gene=self.paths["gene"].stem,
                                            terminator=self.paths["terminator"].stem,
                                            **self.primers["names"])

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


if __name__ == "__main__":

    workdir = "/home/bjorn/Desktop/python_packages/ypkpathway/tests/test/new"

    datafolders = (
        "/home/bjorn/Desktop/python_packages/ypkpathway/tests/test",
        "/home/bjorn/Desktop/YeastPathwayKit/sequences",)

    name = "pTA1_TDH3_ScCTT1_PGI1"

    self = TranscriptionalUnit(name, workdir, datafolders)

    self.find_primers()
    self.copy_primers()
    self.read_primers()

    self.find_sequence_files()
    self.copy_sequence_files()
    self.read_sequence_files()

    self.find_enzymes()

    self.pcr()

    self.format_code()

    self.execute_nb()
