#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pathlib import Path
from pydna.readers import read
from pydna.parsers import parse_primers
# from pydna.amplicon import Amplicon
from pydna.dseqrecord import Dseqrecord
import importlib
import jupytext
import pandas as pd
from Bio.Restriction import ZraI, AjiI, EcoRV
from Bio.Restriction import CommOnly as Co

# https://stackoverflow.com/questions/29850801/subclass-pathlib-path-fails


(p577,
 p578,
 p468,
 p467,
 p567,
 p568,
 p775,
 p778,
 p342) = parse_primers("standard_primers.fasta")


class PathWay:
    """docstring."""

    def __init__(self,
                 pwname: str,
                 kit_folder: Path = None,
                 user_folder: Path = None,
                 work_folder: Path = None,
                 element_backbone_name: str = "pYPKa",
                 tu__backbone_name: str = "pYPKpw"):
        self.pwname = pwname
        self.kit_folder = Path(kit_folder)
        self.user_folder = Path(user_folder)
        self.work_folder = Path(work_folder)
        self.element_backbone_name = element_backbone_name
        self.tu__backbone_name = tu__backbone_name

        backbone, *elements = pwname.split("_")

        if len(elements) % 2 != 1:
            raise ValueError("Uneven number of elements")

        if len(elements) < 3:  # At least one TU
            raise ValueError("At least a promoter, a gene and a terminator")

        self.vectors = [work_folder / Path(pwname).with_suffix(".gb"),
                        kit_folder / Path(backbone).with_suffix(".gb"),
                        kit_folder / Path(element_backbone_name).with_suffix(".gb"),
                        kit_folder / Path(tu__backbone_name).with_suffix(".gb")]
        
        self.TUvectors = []

        for i in range(0, len(elements) - 1, 2):

            p, g, t = elements[i:i+3]

            self.TUvectors.append(
                TranscriptionalUnit(f"{tu__backbone_name}_{p}_{g}_{t}.gb",
                                    pathway=self))
    
    def inspect_project():
        pass
        
    
    
    def generate_project(self):
        pass


def cloned(vector: Dseqrecord, enzymes, candidate: Dseqrecord) -> bool:
    if len(candidate) <= len(vector):
        return False
    candidate2 = str(candidate.seq).lower() * 2
    if enzymes:
        for linear_vector in vector.cut(enzymes):
            if str(linear_vector.seq).lower() in candidate2:
                msg = True
            else:
                msg = False
    else:
        st = candidate2.find("gttcttgtctcattgccacattcataagt")
        en = candidate2.find("gttctgatcctcgagcatcttaagaattc")
        extractbb = candidate2[st:en]
        if extractbb in str(vector.seq).lower()*2:
            msg = True
        else:
            msg = False
    return msg


class VectorPath(Path):

    _flavour = type(Path())._flavour

    def __init__(self,
                 backbone=None,
                 *args,
                 **kwargs):
        super().__init__()
        self.dsr = None
        self.bb = backbone

    def read_Dseqrecord(self):
        try:
            text = self.read_text()
        except FileNotFoundError:
            pass
        else:
            try:
                seq = read(text)
            except ValueError:
                print("can't read")
            else:
                self.dsr = seq
            return bool(self.dsr)

    def check(self):
        elements = self.stem.split("_")
        if len(elements) == 3:  # pYPKa vector
            bb, enz, gene = elements
            self.enz = {"Z": ZraI, "A": AjiI, "E": EcoRV}[enz]
            bb = read(str(bb))
            msg = cloned(bb, self.enz, seq)
        elif len(elements) == 5:  # tu vector
            bb, p, g, t = elements
            if bb == "pYPK0":
                bb = "pYPKpw"
            bb = read(str(bb))
            self.enz = None
            msg = cloned(bb, enz, seq)
        elif len(elements) == 1:  # backbone vector
            mcs = pcr(p577, p578, self.dsr)[len(p577):-len(p578)]
            eb = mcs.single_cutters(Co) & self.dsr.single_cutters(Co)
            _, self.enz,  = sorted([(abs(len(f1)-len(f2)), e) for f1, f2 in [mcs.cut(e) for e in eb]])
        else:
            raise ValueError("name has wrong format.")




class TranscriptionalUnit:
    """docstring.

    Single-Gene Expression Vector
    """

    def __init__(self,
                 name: str = "",
                 pathway: PathWay = None):

        self.name = name
        self.pathway = pathway
        self.path = None

        (self.backbone_vector_name,
         self.promoter_vector_name,
         self.gene_vector_name,
         self.terminator_vector_name) = name.removesuffix(".gb").split("_")

    def find_tu_file(self):
        """docstring."""
        uf, kf, wf = [self.pathway.user_folder,
                      self.pathway.kit_folder,
                      self.pathway.work_folder]
        for f in wf, uf, kf, wf:
            self.path = f/Path(self.name)
            if self.path.exists():
                break
        return Path(self.name).exists()

    def inspect_tu_file(self):
        """docstring."""
        result = cloned()
        return result

    def find_te_files(self):
        """docstring."""
        kf = self.pathway.kit_folder
        uf = self.pathway.user_folder
        wf = self.pathway.work_folder

        self.backbone_vector_exists = False
        self.promoter_vector_exists = False
        self.gene_vector_exists = False
        self.terminator_vector_exists = False

        for f in [kf, uf, wf, kf]:
            self.backbone_vector_path = f/VectorPath(
                f"{self.backbone_vector_name}").with_suffix(".gb")
            if self.backbone_vector_path.exists():
                self.backbone_vector_exists = True
                break

        for f in [kf, uf, wf, kf]:
            self.promoter_vector_path = f/VectorPath(
                f"pYPKa_Z_{self.promoter_vector_name}").with_suffix(".gb")
            if self.promoter_vector_path.exists():
                self.promoter_vector_exists = True
                break

        for f in [wf, uf, kf, wf]:
            self.gene_vector_path = f/VectorPath(
                f"pYPKa_A_{self.gene_vector_name}").with_suffix(".gb")
            if self.gene_vector_path.exists():
                self.gene_vector_exists = True
                break

        for f in [kf, uf, wf, kf]:
            self.terminator_vector_path = f/VectorPath(
                f"pYPKa_E_{self.terminator_vector_name}.gb").with_suffix(".gb")
            if self.terminator_vector_path.exists():
                self.terminator_vector_exists = True
                break















def df(name):
    f = importlib.resources.files("ypkpathway") / Path("data") / name
    return f.read_text()


def load_files(pth="GeneticAssembly"):
    pth = Path(pth)
    pth.mkdir(parents=True, exist_ok=True)
    return list(pth.glob("*.gb"))


def process_pathway_name(pwname: str,
                         kit_folder: Path,
                         user_folder: Path,
                         work_folder: Path,
                         el_bb_name: str = "pYPKa",
                         tu_bb_name: str = "pYPKpw") -> list:
    """
    Process pathway name.

    A list of Path objects for the three folders.
    """
    backbone, *elements = pwname.split("_")

    assert len(elements) % 2 == 1  # Uneven number of elements

    assert len(elements) >= 3  # At least one TU

    kit_folder = Path(kit_folder)
    user_folder = Path(user_folder)
    work_folder = Path(work_folder)
    
    # tu_bb_name = "pYPK0" if tu_bb_name=="pYPKpw" else tu_bb_name

    vectors = [
        work_folder / Path(pwname).with_suffix(".gb"),
        kit_folder / Path(backbone).with_suffix(".gb"),
        kit_folder / Path(el_bb_name).with_suffix(".gb"),
        kit_folder / Path(tu_bb_name).with_suffix(".gb")]

    for i in range(0, len(elements) - 1, 2):

        p, g, t = elements[i:i+3]

        vectors.append(user_folder / Path(f"{tu_bb_name}_{p}_{g}_{t}.gb"))
        vectors.append(kit_folder / Path(f"pYPKa_Z_{p}.gb"))
        vectors.append(user_folder / Path(f"pYPKa_A_{g}.gb"))
        vectors.append(kit_folder / Path(f"pYPKa_E_{t}.gb"))

    return vectors


def check_paths(vectors: list[Path]) -> df:
    """docstring."""
    exists = []
    readable = []
    length = []
    circular = []
    clonetype = []

    vectordict = {v.stem: v for v in vectors}

    for v in vectors:
        if v.exists():
            exists.append(True)
            try:
                seq = read(v)
            except ValueError:
                readable.append(False)
                length.append(None)
                circular.append(None)
                clonetype.append(None)
            else:
                readable.append(True)
                if len(seq) < 3128:
                    length.append(False)
                    circular.append(None)
                    clonetype.append(None)
                else:
                    length.append(True)
                    if seq.linear:
                        circular.append(False)
                        clonetype.append(None)
                    else:
                        circular.append(True)
                        elements = v.stem.split("_")
                        if len(elements) == 3:  # pYPKa vector
                            bb, enz, gene = elements
                            enz = {"Z": ZraI, "A": AjiI, "E": EcoRV}[enz]
                            bb = vectordict[bb]
                            bb = read(bb)
                            clonetype.append(cloned(bb, enz, seq))
                        elif len(elements) == 5:  # tu vector
                            bb, p, g, t = elements
                            if bb == "pYPK0": 
                                bb = "pYPKpw"
                            bb = vectordict[bb]
                            bb = read(bb)
                            enz = None
                            clonetype.append(cloned(bb, enz, seq))
                        elif len(elements) == 1:  # backbone vector
                            clonetype.append(None)             
        else:
            exists.append(False)
            readable.append(None)
            length.append(None)
            circular.append(None)
            clonetype.append(None)
            
    # breakpoint()
    df = pd.DataFrame(zip(*[vectors, 
                            exists, 
                            readable, 
                            length, 
                            circular, 
                            clonetype]))
    df.columns = "Path exists readable length circular clonetype".split()

    return df

def assemble_tu_vectors(tu_vector_names, usergbfiles, kitfgbfiles):
    files = {}
    for tu_name in tu_vector_names:
        if tu_name.name in [f.name for f in usergbfiles]:
            files[tu_name] = tu_name.read_text()
        else:
            files[tu_name] = TU_assembly(tu_name)
    return files


def TU_assembly(tu_name):

    tu_name = Path(tu_name)

    files = {}

    bb_name, promoter_name, gene_name, terminator_name = tu_name.stem.split(
        "_"
    )

    promoter_vector_path = kit_folder / Path(f"pYPKa_Z_{promoter_name}.gb")
    gene_vector_path = user_folder / Path(f"pYPKa_A_{gene_name}.gb")
    terminator_vector_path = kit_folder / Path(f"pYPKa_E_{terminator_name}.gb")

    mdnbtemplate = df("nb_template_pYPK0_tp_gene_tp.md")

    formattedmd = mdnbtemplate.format(
        tpz=promoter_name, gene=gene_name, tpe=terminator_name
    )

    ipynb = jupytext.reads(formattedmd, fmt="md")

    files[tu_name.with_suffix(".md")] = jupytext.writes(ipynb, fmt="ipynb")

    return files


def execute_notebooks(files, log, working_dir="ypkassembly", print=print):

    try:
        Path(working_dir).mkdir(parents=True)
    except OSError as exception:
        if exception.errno == errno.EEXIST:
            print(
                "The {} directory already exists! Please delete or choose another name.".format(
                    working_dir
                )
            )
        else:
            print("The {} directory could not be created".format(working_dir))
        return None, None

    cwd = os.getcwd()

    msg = "created subdirectory {}\n".format(working_dir)

    print(msg)

    log += msg

    os.chdir(working_dir)

    msg = "saving files sequence files and images.\n"
    print(msg)
    log += msg

    notebooks = []
    obj = notedown.MarkdownReader()

    for name, content in files.items():
        if name.endswith(".md"):
            newname = os.path.splitext(name)[0] + ".ipynb"
            msg = "\nsaving: " + newname
            print(msg)
            log += msg
            nb = nbformat.write(obj.to_notebook(content), newname)
            notebooks.append(newname)
        else:
            mode = {True: "wb", False: "w"}[hasattr(content, "decode")]
            with open(name, mode) as f:
                f.write(content)

    del files

    ep = ExecutePreprocessor()
    ep.timeout = 1200  # seconds
    ep.interrupt_on_timeout = True

    print("\n")
    log += "\n"

    msg = "\nexecuting pYPKa notebooks..\n"
    print(msg)
    log += msg

    shell = InteractiveShell.instance()

    g = {}
    l = {}
    resources = {}

    for name in notebooks[:-1]:
        msg = "\nexecuting: " + name
        print(msg)
        log += msg
        with open(name, "r", encoding="utf-8") as f:
            nb = nbformat.read(f, 4)
        nb_executed, resources = ep.preprocess(nb, resources=resources)
        for cell in nb.cells:
            if cell.cell_type == "code":
                code = shell.input_transformer_manager.transform_cell(
                    cell.source
                )
                exec(code, g, l)
        nbformat.write(nb, name)
        g = {}
        l = {}

    os.chdir(cwd)

    fl = FileLink(os.path.join(working_dir, "pw.ipynb"))

    return fl, log


if __name__ == "__main__":

    # name = "pTA1_TDH3_ScATF1_PGI1"
    # name = "pTA1_PDC1_EcfabH_TEF1_EcfabD_FBA1_EcfabG_RPL22A_EcacpP_TDH3_EcfabF_UTR2_EcfabB_TPI1_EcfabA_PMP3_EcfabZ_ENO2_Athmod1_RPL5_AthfatA1_RPL16A_AthfatB_RPL17A_EcacpH_RPL16B_EcacpS_TMA19"

    pw = PathWay("pTA1_TDH3_ScATF1_PGI1_ScCTT1_TEF1",  
                 kit_folder = Path("YeastPathwayKit/sequences"),
                 user_folder = Path("UserSequences"),
                 work_folder = Path("work"))

    
    x, y = pw.TUvectors
    
    x.find_tu_file()


    












