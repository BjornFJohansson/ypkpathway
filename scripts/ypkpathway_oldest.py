#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ypkpathway_cli

Usage: ypkpathway from_file <path> [<dir>] [--no_pYPKa_A]
       ypkpathway from_name <name> [<dir>]
       ypkpathway -h|--help
       ypkpathway -v|--version

Arguments:
    <path>  path to data file containing sequences to be assembled

    <dir>   Optional directory to put generated sequence files, defaults to
            <ypk_assembly> in the current working directory.

Options:
    -h, --help      Show this screen.
    -v, --version   Show version.
"""


import io
import re
import sys
import os
import errno
import shutil
import docopt
import subprocess
import itertools
import appdirs
import pathlib

from pkg_resources import resource_filename

import pydna
from pydna.readers import read
from pydna.parsers import parse


import nbformat

from nbconvert.preprocessors.execute import ExecutePreprocessor

from IPython.core.interactiveshell import InteractiveShell

from IPython.display import FileLink, FileLinks

import notedown

re_cas  = re.compile(r"pYPK0_([^\d\W]\w{2,15})_([^\d\W]\w{2,15})_([^\d\W]\w{2,15})")
re_cas  = re.compile(r"pYPK0_([^_]{2,15})_([^_]{2,15})_([^_]{2,15})")
re_Z    = re.compile(r"pYPKa_Z_([^\d\W]\w{2,15})")
re_A    = re.compile(r"pYPKa_A_([^\d\W]\w{2,15})")
re_E    = re.compile(r"pYPKa_E_([^\d\W]\w{2,15})")

def cloned(vector, enzyme, candidate):
    if len(candidate) <= len(vector):
        return 0
    candidate2 = str(candidate.seq).lower()*2
    for linear_vector in vector.cut(enzyme):
        if str(linear_vector.seq).lower() in candidate2:
            return len(candidate) - len(vector)
    return 0

def read_data_file(name):
    with open( resource_filename("ypkpathway", os.path.join("data", name)), "r") as f:
        data = f.read()
    return data


def read_bin_file(name):
    with open( resource_filename("ypkpathway", os.path.join("data", name)), "rb") as f:
        data = f.read()
    return data



from typing import Any

def pathway(text    : Any,
            dir_    : str = "ypkassembly",
            YPKa_A  : bool = True,
            print = print):

    nb, log = None, None

    sequence_list = parse(text)

    if text.startswith("pYPK0_") and not sequence_list:
        nb, log = from_name()
    elif sequence_list:
        nb, log = from_list(sequence_list, dir_, pYPKa_A = pYPKa_A)

    return nb, log

def from_name(pwname):

    files = {"standard_primers.txt"     : read_data_file("standard_primers.txt"),
             "pYPKa.gb"                 : read_data_file("pYPKa.gb"),
             "pYPKpw.gb"                : read_data_file("pYPKpw.gb"),
             "tp_g_tp.png"              : read_bin_file("tp_g_tp.png"),
             "pYPK_ZE.png"              : read_bin_file("pYPK_ZE.png"),
             "pYPK_A.png"               : read_bin_file("pYPK_A.png"),
             "pw.png"                   : read_bin_file("pw.png"),
             "start.bat"                : read_data_file("start.bat"),
             "start.sh"                 : read_data_file("start.sh"),}

    cas_vectors = []
    tp_gene_tp_links = ""
    pYPKa_clones=""
    log=""

    log+= f"filename = {pwname}\n\n"

    items = [n.strip().split("tp",1)[0] for n in pwname.split("_")[1:] ]

    log+= "\n".join(items)

    promoters, genes, terminators = items[0:-1:2],items[1::2],items[2::2]

    log+= "\npromoters"
    log+=str(promoters)
    log+= "\ngenes"
    log+=str(genes)
    log+= "\nterminators"
    log+=str(terminators)

    appdir = pathlib.Path(appdirs.user_data_dir("ypkpathway")).joinpath("constructs")

    log+=f"appdir {appdir}"

    cwd = pathlib.Path().cwd()

    log+=f"cwd {cwd}"

    existing_files  = { u.name:u for u in itertools.chain(appdir.rglob("pYPK*.gb"), cwd.rglob("pYPK*.gb")) }

    log+=f"existing_files {len(existing_files)}"

    missing = []

    for p,g,t in zip(promoters, genes, terminators):

        cassette_vector_path = pathlib.Path( f"pYPK0_{p}_{g}_{t}.gb" )

        cas_vectors.append( cassette_vector_path.name )

        tp_gene_tp_links += "\n[{}]({})\n".format( cassette_vector_path.name, cassette_vector_path )

        if cassette_vector_path.name in existing_files:
            files[cassette_vector_path.name] = existing_files[cassette_vector_path.name].read_text()
            continue
        else:
            promoter_vector_path   = f"pYPKa_Z_{p}.gb"
            gene_vector_path       = f"pYPKa_A_{g}.gb"
            terminator_vector_path = f"pYPKa_E_{t}.gb"

            for item in promoter_vector_path, gene_vector_path, terminator_vector_path:
                if item in existing_files:
                    files[item] = existing_files[item].read_text()
                else:
                    missing.append(item)

            nbtemp = read_data_file("nb_template_pYPK0_tp_gene_tp.md")

            files[cassette_vector_path.with_suffix(".md").name] = nbtemp.format(tpz=p, gene=g,tpe=t)

    nbtemp = read_data_file("nb_template_pYPK0_pw_from_name.md")

    cas_vectors_code = "".join(f'"{v}",' for v in cas_vectors)


    files["pw.md"] = nbtemp.format( name=pwname,
                                    pwname=pwname,
                                    tp_gene_tp_links = tp_gene_tp_links,
                                    cas_vectors=cas_vectors_code,
                                    length=len(genes))

    if missing:
        log += "\n\n".join(str(missing))
    else:
        log += "\nOK."

    return files, log



def from_list(sequence_list   : list,
                      pYPKa_A : bool = True):


    try:
        with open(pth, "r") as f:
            text=f.read()
    except IOError:
        print(pth, "could not be opened!")
        sys.exit(1)

    seqs = parse(text)

    if len(seqs)==0: # text has to contain some sequences
        print("No sequences in argument.")
        return None, None

    log=""

    pYPK0 = read(read_data_file("pYPK0.gb"))
    pYPKa = read(read_data_file("pYPKa.gb"))

    from Bio.Restriction import ZraI, AjiI, EcoRV

    files = {"standard_primers.txt"     : read_data_file("standard_primers.txt"),
             "pYPKa.gb"                 : read_data_file("pYPKa.gb"),
             "pYPKpw.gb"                : read_data_file("pYPKpw.gb"),
             "tp_g_tp.png"              : read_bin_file("tp_g_tp.png"),
             "pYPK_ZE.png"              : read_bin_file("pYPK_ZE.png"),
             "pYPK_A.png"               : read_bin_file("pYPK_A.png"),
             "pw.png"                   : read_bin_file("pw.png"),
             "start.bat"                : read_data_file("start.bat"),
             "start.sh"                 : read_data_file("start.sh"),}

    files[f"INDATA_{pth.stem}.txt"] = text

    cas_vectors = ""
    tp_gene_tp_links = ""
    pYPKa_clones=""
    pwname = "pYPK0"
    genes = 0
    nbflag=False

    while pth:
        genes+=1
        first = pth.pop(0)
        # is sequence a tp-gene-tp vector?
        if cloned(pYPK0, (ZraI, EcoRV),  first):
            m = re_cas.search(first.description)
            if not m:
                raise Exception( "{} is a pYPK0 tp-gene_tp sequence but it was not correctly named.".format(first.description))
            fn = first.description+".gb"
            files[fn] = first.format("gb")
            cas_vectors+= fn+"\n"
            tp_gene_tp_links+= "\n[{}]({})\n".format( first.description, fn )
            tp1_description  = m.group(1)
            gene_description = m.group(2)
            tp2_description  = m.group(3)
            genes+=1
        else:
            try:
                middle = pth.pop(0)
                last   = pth.pop(0)
            except IndexError:
                raise Exception("not enough sequences")

            prom, gene, term = first, middle, last

            if cloned(pYPKa, ZraI,  prom):
                m = re_Z.search(prom.description)
                if not m:
                    raise Exception( "{} is a pYPKa_A_gene sequence but was incorrectly named.".format(gene.description))
                prom_description = m.group(1)
                files[m.group(0)+".gb"] = prom.format("gb")
            else:
                #print("Z"+str(files.has_key("pYPKa_ZE_{}.md".format(prom.id)))+prom.id)
                if "pYPKa_ZE_{}.md".format(prom.id) not in files:
                    files[prom.id+".gb"] = prom.format("gb")
                    nbtemp = read_data_file("nb_template_pYPKa_ZE_insert.md")
                    files["pYPKa_ZE_{}.md".format(prom.id)] = nbtemp.format(tp=prom.id)
                    pYPKa_clones+="[pYPKa_ZE_{n}](pYPKa_ZE_{n}.ipynb)  \n".format(n=prom.id)
                prom_description = prom.id

            if cloned(pYPKa, AjiI,  gene):
                m = re_A.search(gene.description)
                if not m:
                    raise Exception( "{} is a pYPKa_A_gene sequence but was incorrectly named.".format(gene.description))
                gene_description = m.group(1)
                files[m.group(0)+".gb"] = gene.format("gb")
                if not pYPKa_A:
                    nbflag=True

            else:
                n = "pYPKa_A_{}".format(gene.locus)
                files[gene.locus+".gb"] = gene.format("gb")
                if pYPKa_A:
                    nbtemp = read_data_file("nb_template_pYPKa_A_insert.md")
                    files[n+".md"] = nbtemp.format(insert=gene.locus)
                    gene_description = gene.locus
                    pYPKa_clones+="[{}]({}.ipynb)  \n".format(n, n)
                else:
                    gene_description = gene.locus

            if cloned(pYPKa, EcoRV, term):
                m = re_E.search(term.description)
                if not m:
                    raise Exception( "{} is a pYPKa_A_gene sequence but was incorrectly named.".format(gene.description))
                term_description = m.group(1)
                files[m.group(0)+".gb"] = term.format("gb")
            else:
                #print("E"+str(files.has_key("pYPKa_ZE_{}.md".format(term.id)))+term.id)
                if "pYPKa_ZE_{}.md".format(term.id) not in files:
                    files[term.id+".gb"] = term.format("gb")
                    nbtemp = read_data_file("nb_template_pYPKa_ZE_insert.md")
                    files["pYPKa_ZE_{}.md".format(term.id)] = nbtemp.format(tp=term.id)
                    pYPKa_clones+="[pYPKa_ZE_{n}](pYPKa_ZE_{n}.ipynb)  \n".format(n=term.id)
                term_description = term.id

            x = "pYPK0_{}_{}_{}".format(prom_description, gene_description, term_description)

            if pYPKa_A or nbflag:
                nbtemp = read_data_file("nb_template_pYPK0_tp_gene_tp.md")
                files[x+".md"] = nbtemp.format(tpz=prom_description,
                                                gene=gene_description,
                                                tpe=term_description)
            else:
                nbtemp = read_data_file("nb_template_pYPK0_tp_gene_tp_gap_repair.md")
                files[x+".md"] = nbtemp.format(tpz=prom_description,
                                                gene=gene.locus,
                                                tpe=term_description)
            nbflag=False

            cas_vectors+="\n"+x+".gb\n"
            tp_gene_tp_links+="[{}]({}.ipynb)  \n".format(x, x)

            pwname+="_{}".format(gene_description)

        nbtemp = read_data_file("nb_template_pYPK0_pw_from_sequence.md")

        files["pw.ipynb"] = nbtemp.format(name=pwname,
                                          filename=os.path.basename(dir_),
                                          tp_gene_tp_links = tp_gene_tp_links,
                                          cas_vectors=add_space(cas_vectors, 17),
                                          pYPKa_clones=pYPKa_clones,
                                          length=len(genes))

    return files, log

    ###########################################################################



def execute_notebooks(files,log, working_dir="ypkassembly", print=print):

    try:
        pathlib.Path(working_dir).mkdir(parents=True)
    except OSError as exception:
        if exception.errno == errno.EEXIST:
            print("The {} directory already exists! Please delete or choose another name.".format(working_dir))
        else:
            print("The {} directory could not be created".format(working_dir))
        return None, None


    cwd = os.getcwd()

    msg = "created subdirectory {}\n".format(working_dir)

    print(msg)

    log+=msg

    os.chdir(working_dir)

    msg = "saving files sequence files and images.\n"
    print(msg)
    log+=msg

    notebooks = []
    obj = notedown.MarkdownReader()

    for name, content in files.items():
        if name.endswith(".md"):
            newname = os.path.splitext(name)[0]+".ipynb"
            msg = "\nsaving: "+newname
            print(msg)
            log+=msg
            nb = nbformat.write(obj.to_notebook(content), newname)
            notebooks.append(newname)
        else:
            mode = {True:"wb", False:"w"}[hasattr(content, "decode")]
            with open(name, mode) as f:
                f.write(content)

    del files

    ep = ExecutePreprocessor()
    ep.timeout = 1200              # seconds
    ep.interrupt_on_timeout = True

    print("\n")
    log+="\n"

    msg = "\nexecuting pYPKa notebooks..\n"
    print(msg)
    log+=msg

    shell = InteractiveShell.instance()

    g={}
    l={}
    resources = {}
    
    for name in notebooks[:-1]:
        msg = "\nexecuting: "+name
        print(msg)
        log+=msg
        with open(name, 'r', encoding='utf-8') as f:
            nb = nbformat.read(f, 4)
        nb_executed, resources = ep.preprocess(nb,
                                               resources=resources)
        for cell in nb.cells:
            if cell.cell_type == 'code':
                code = shell.input_transformer_manager.transform_cell(cell.source)
                exec(code, g, l)
        nbformat.write(nb, name)
        g={}
        l={}

    os.chdir(cwd)

    fl = FileLink(os.path.join(working_dir, "pw.ipynb"))

    return fl, log

def main():

    import appdirs
    import pathlib

    appdir = pathlib.Path(appdirs.user_data_dir("ypkpathway")).joinpath("constructs")

    appdir.mkdir(parents=True, exist_ok=True)

    try:
        arguments = docopt.docopt(__doc__)
    except docopt.DocoptExit as e:
        sys.exit(1)

    working_dir = "ypk_assembly"

    if arguments["<dir>"]:
        working_dir = str(arguments["<dir>"])

    if arguments["--no_pYPKa_A"]:
        pYPKa_A = False
    else:
        pYPKa_A = True

    if arguments["--version"]:
        from ._version import get_versions
        __version__ = get_versions()["version"][:5]
        del get_versions
        print("ypkpathway version:",__version__)
        print("     pydna version:",pydna.__version__)


    if arguments["<path>"]:

        pth = pathlib.Path( str(arguments["<path>"]) )

        print("Assembly started. (This might take a while...)")

        files, log = from_list( pth, pYPKa_A=pYPKa_A )

    elif arguments["<name>"]:

        print("Assembly started. (This might take a while...)")

        files, log = from_name( str(arguments["<name>"]))

    if log.endswith("OK."):

        execute_notebooks(files,log, working_dir)

        print("opening IPython notebook")

        subprocess.Popen(["jupyter", "notebook", os.path.join(working_dir, "pw.ipynb")])
    else:
        print(log)


def pYPKa_ZE_ipynb_generator(tp, dir_="pYPKa_ZE_vectors"):

    cwd = os.getcwd()

    try:
        os.makedirs(dir_)
    except OSError as exception:
        if exception.errno == errno.EEXIST:
            passdir_
        else:
            print("The {} directory could not be created".format(dir_))
            return None

    os.chdir(dir_)

    with open("standard_primers.txt","w") as f: f.write(read_data_file("standard_primers.txt"))
    with open("pYPKa.gb","w") as f: f.write(read_data_file("pYPKa.gb"))
    with open("pYPK_ZE.png","w") as f: f.write(read_bin_file("pYPK_ZE.png"))
    with open(tp.id+".gb","w") as f: f.write(tp.format("gb"))

    #nbtemp = execute_notebooks(files,log, working_dir)d_data_file("nb_template_pYPKa_ZE_insert.md"))

    name = "pYPKa_ZE_{}.ipynb".format(tp.id)

    obj = notedown.MarkdownReader()

    nb = obj.to_notebook(nbtemp.format(tp=tp.id))

    pp = ExecutePreprocessor()
    pp.timeout = 120 # seconds
    pp.interrupt_on_timeout = True

    shell = InteractiveShell.instance()

    nb_executed, resources = pp.preprocess(nb, resources={})

    g={}
    l={}

    old_stdout = sys.stdout
    redirected_output = sys.stdout = io.StringIO()

    for cell in nb.cells:
        if cell.cell_type == 'code':
            code = shell.input_transformer_manager.transform_cell(cell.source)
            exec(code, g, l)

    sys.stdout = old_stdout

    nbformat.write(nb, name)

    os.chdir(cwd)

    return FileLinks(dir_)

if __name__ == "__main__":
    main()


