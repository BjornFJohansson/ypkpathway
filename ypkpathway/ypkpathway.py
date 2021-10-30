#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ypkpathway

Usage: ypkpathway assembly_from_file  <path> [<dir>] [--no_pYPKa_A]
       ypkpathway assembly_from_name  <name> [<dir>]
       ypkpathway pYPKa_A_from_csv    <path> [<dir>] --email=<email>
       ypkpathway pYPK0_A_from_csv    <path> [<dir>] --email=<email>
       ypkpathway pYPKa_ZE_from_csv   <path> [<dir>] --email=<email>
       ypkpathway -h|--help
       ypkpathway -v|--version

Arguments:
    <path>  path to data file containing sequences to be assembled or
            data in specific csv format.

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
import docopt
#import itertools
import appdirs
import pathlib
import csv

from pkg_resources import resource_filename
from datetime import datetime
from glob import glob
from typing import Any

from tqdm import tqdm
import pydna
from pydna.readers import read
from pydna.parsers import parse
import notedown
import nbformat
from nbconvert.preprocessors.execute import ExecutePreprocessor
from Bio.Restriction import ZraI, EcoRV, AjiI

re_cas  = re.compile(r"pYPK0_([^\d\W]\w{2,15})_([^\d\W]\w{2,15})_([^\d\W]\w{2,15})")
re_cas  = re.compile(r"pYPK0_([^_]{2,15})_([^_]{2,15})_([^_]{2,15})")
re_Z    = re.compile(r"pYPKa_Z_([^\d\W]\w{2,15})")
re_A    = re.compile(r"pYPKa_A_([^\d\W]\w{2,15})")
re_E    = re.compile(r"pYPKa_E_([^\d\W]\w{2,15})")


def print2(text):
    print(text)
    return text


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


def generate(text           : Any,
             pYPKa_A        : bool = True,
             print = print):

    files, error = None, None

    import pydna
    from datetime import datetime
    print(str(datetime.now()))
    import socket
    print(f"host name: {socket.gethostname():20}")
    import getpass
    print(f"user name: {getpass.getuser():20}")
    print(f"pydna version: {pydna.__version__}")
    from IPython import __version__ as IPython_version
    print(f"IPython version: {IPython_version}")
    #jupyter lab --version
    #jupyter notebook --version

    sequence_list = parse(text)

    if text.startswith("pYPK0_") and not sequence_list:
        print("argument is a name")
        files, error = from_name(text, print = print)
    elif sequence_list:
        print(f"argument is a set of {len(sequence_list)} sequences")
        files, error = from_list(sequence_list,
                                  pYPKa_A = pYPKa_A,
                                  print = print)

    return files, error


def from_name(pwname, print = print):

    files = {"standard_primers.txt"     : read_data_file("standard_primers.txt"),
             "pYPKa.gb"                 : read_data_file("pYPKa.gb"),
             "pYPKpw.gb"                : read_data_file("pYPKpw.gb"),
             "tp_g_tp.png"              : read_bin_file("tp_g_tp.png"),
             "pYPK_ZE.png"              : read_bin_file("pYPK_ZE.png"),
             "pYPK_A.png"               : read_bin_file("pYPK_A.png"),
             "pw.png"                   : read_bin_file("pw.png"),
             "start.bat"                : read_data_file("start.bat"),
             "start.sh"                 : read_data_file("start.sh"),}

    error = False
    cas_vectors = []
    tp_gene_tp_links = ""

    print("from_name")

    print(f"filename = {pwname}")

    items = [n.strip().split("tp",1)[0] for n in pwname.split("_")[1:] ]

    promoters, genes, terminators = items[0:-1:2],items[1::2],items[2::2]

    print("promoters")
    print(promoters)
    print("genes")
    print(genes)
    print("terminators")
    print(terminators)

    appdir = pathlib.Path(appdirs.user_data_dir("ypkpathway")).joinpath("sequences")

    print(f"appdir {appdir}")

    cwd = pathlib.Path().cwd()

    print(f"cwd {cwd}")

    existing_files  = { u.name:u for u in (pathlib.Path(p) for p in glob( str(appdir.joinpath("**/pYPK*.gb")), recursive = True)) }

    for p in cwd.glob("pYPK*.gb"):
        existing_files[p.name] = p

    print(f"existing_files {len(existing_files)}")

    obj = notedown.MarkdownReader()

    for p,g,t in zip(promoters, genes, terminators):

        cassette_vector_path = pathlib.Path( f"pYPK0_{p}_{g}_{t}.gb" )

        print(cassette_vector_path.name)

        cas_vectors.append( cassette_vector_path.name )

        tp_gene_tp_links += "\n[{}]({})\n".format( cassette_vector_path.with_suffix(".ipynb").name, cassette_vector_path.with_suffix(".ipynb") )

        if cassette_vector_path.name in existing_files:
            files[cassette_vector_path.name] = existing_files[cassette_vector_path.name].read_text()
            print("found")
            continue
        else:
            promoter_vector_paths   =  f"pYPKa_Z_{p}.gb", f"pYPK0_Z_{p}.gb"
            f=None
            for promoter_vector in promoter_vector_paths:
                try:
                    f = existing_files[promoter_vector].read_text()
                except KeyError:
                    print(f"promoter vector {promoter_vector} not found")
                else:
                    print(f"promoter vector {promoter_vector} found")
                    files[promoter_vector] = f
                    break
            if not f: error = True
            gene_vector_paths =  f"pYPKa_A_{g}.gb", f"pYPK0_A_{g}.gb"
            f=None
            for gene_vector in gene_vector_paths:
                try:
                    f = existing_files[gene_vector].read_text()
                except KeyError:
                    print(f"gene vector {gene_vector} not found")
                else:
                    print(f"gene vector {gene_vector} found")
                    files[gene_vector] = f
                    break
            if not f: error = True
            terminator_vector_paths =  f"pYPKa_E_{t}.gb", f"pYPK0_E_{t}.gb"
            f=None
            for terminator_vector in terminator_vector_paths:
                try:
                    f = existing_files[terminator_vector].read_text()
                except KeyError:
                    print(f"terminator vector {terminator_vector} not found")
                else:
                    print(f"terminator vector {terminator_vector} found")
                    files[terminator_vector] = f
                    break
            if not f: error = True

            nbtemp = read_data_file("nb_template_pYPK0_tp_gene_tp.md")

            files[cassette_vector_path.with_suffix(".ipynb").name] = obj.to_notebook(nbtemp.format(tpz=p,
                                                                                                   gene=g,
                                                                                                   tpe=t,
                                                                                                   promoter_vector=promoter_vector,
                                                                                                   gene_vector=gene_vector,
                                                                                                   terminator_vector=terminator_vector))

    nbtemp = read_data_file("nb_template_pYPK0_pw_from_name.md")

    cas_vectors_code = "".join(f'"{v}",' for v in cas_vectors)

    files["pw.ipynb"] = obj.to_notebook( nbtemp.format( name=pwname,
                                         pwname=pwname,
                                         tp_gene_tp_links = tp_gene_tp_links,
                                         cas_vectors=cas_vectors_code,
                                         length=len(genes)) )

    return files, error


def from_list(seqs        : list,
              pYPKa_A     : bool = True,
              print = print):

    if len(seqs)==0:
        print("No sequences in argument.")
        return {}, True

    obj = notedown.MarkdownReader()

    pYPK0 = read(read_data_file("pYPK0.gb"))
    pYPKa = read(read_data_file("pYPKa.gb"))

    files = {"standard_primers.txt"     : read_data_file("standard_primers.txt"),
             "pYPKa.gb"                 : read_data_file("pYPKa.gb"),
             "pYPKpw.gb"                : read_data_file("pYPKpw.gb"),
             "tp_g_tp.png"              : read_bin_file("tp_g_tp.png"),
             "pYPK_ZE.png"              : read_bin_file("pYPK_ZE.png"),
             "pYPK_A.png"               : read_bin_file("pYPK_A.png"),
             "pw.png"                   : read_bin_file("pw.png"),
             "start.bat"                : read_data_file("start.bat"),
             "start.sh"                 : read_data_file("start.sh"),}

    error = False

    files["INDATA.txt"] = "\n\n".join(s.format() for s in seqs)

    cas_vectors = []
    tp_gene_tp_links = ""
    pYPKa_clones=""

    nbflag=False
    promoters,genes,terminators = [],[],[]

    while seqs:
        first = seqs.pop(0)
        if cloned(pYPK0, (ZraI, EcoRV),  first):
            print(f"{first.description} sequence is a promoter-gene-terminator vector.")
            m = re_cas.search(first.description)
            if not m:
                print(f"{first.description} is a pYPK0 tp-gene_tp sequence but it was not correctly named.")
                error = True
            fn = f"{first.description}.gb"
            files[fn] = first.format("gb")
            cas_vectors.append( fn )
            print(f"{fn} added to file list")
            #tp_gene_tp_links+= "\n[{}]({})\n".format( first.description, fn )
            promoters.append( m.group(1)[:-2] if m.group(1).endswith("tp") else m.group(1))
            genes.append( m.group(2)[:-2] if m.group(2).endswith("tp") else m.group(2))
            terminators.append( m.group(3)[:-2] if m.group(3).endswith("tp") else m.group(3))
            print(f"promoters {' '.join(promoters)}")
            print(f"genes {' '.join(genes)}")
            print(f"terminators {' '.join(terminators)}")
        else:
            try:
                middle = seqs.pop(0)
                last   = seqs.pop(0)
            except IndexError:
                print("not enough sequences")
                error = True

            prom, gene, term = first, middle, last

            if cloned(pYPKa, ZraI,  prom):
                print(f"{prom.description} sequence is a pYPKa_Z_promoter sequence.")
                m = re_Z.search(prom.description)
                if not m:
                    print(f"{prom.description} was incorrectly named.")
                    print("this has to be fixed before assembly can proceed.")
                    error = True
                prom_description = m.group(1)[:-2] if m.group(1).endswith("tp") else m.group(1)
                files[f"pYPKa_Z_{prom_description}.gb"] = prom.format("gb")
                print(f"pYPKa_Z_{prom_description}.gb added to file list.")
            else:
                print(f"{prom.description} is a linear DNA fragment.")
                prom_description = prom.id[:-2] if prom.id.endswith("tp") else prom.id
                tp_vector = f"pYPKa_ZE_{prom_description}.ipynb"
                if tp_vector not in files:
                    print(f"{tp_vector} not found.")
                    print(f"{tp_vector} will be designed")
                    files[prom_description+".gb"] = prom.format("gb")
                    nbtemp = read_data_file("nb_template_pYPKa_ZE_insert.md")
                    files[tp_vector] = obj.to_notebook( nbtemp.format(tp=prom_description) )
                    pYPKa_clones+="[pYPKa_ZE_{n}](pYPKa_ZE_{n}.ipynb)\n".format(n=prom_description)
                    print(f"{tp_vector} added to file list.")
                else:
                    print(f"{tp_vector} already available.")

            if cloned(pYPKa, AjiI,  gene):
                print(f"{gene.description} sequence is a pYPKa_A_gene sequence.")
                m = re_A.search(gene.description)
                if not m:
                    print(f"{gene.description} was incorrectly named.")
                    print("this has to be fixed before assembly can proceed.")
                    error = True
                gene_description = m.group(1)
                files[f"pYPKa_A_{gene_description}.gb"] = gene.format("gb")
                print(f"pYPKa_A_{gene_description}.gb added to file list.")
                if not pYPKa_A: nbflag=True
            else:
                print(f"{gene.description} is a linear DNA fragment.")
                n = "pYPKa_A_{}".format(gene.locus)
                files[gene.locus+".gb"] = gene.format("gb")
                if pYPKa_A:
                    nbtemp = read_data_file("nb_template_pYPKa_A_insert.md")
                    files[n+".ipynb"] = obj.to_notebook( nbtemp.format(insert=gene.locus) )
                    print(f"{n}.ipynb added to file list.")
                    gene_description = gene.locus
                    pYPKa_clones+="[{}]({}.ipynb)\n".format(n, n)
                else:
                    print(f"{gene.locus} will be cloned by gap repair.")
                    gene_description = gene.locus

            if cloned(pYPKa, EcoRV, term):
                print(f"{term.description} sequence is a pYPKa_E_terminator sequence.")
                m = re_E.search(term.description)
                if not m:
                    print(f"{gene.description} was incorrectly named.")
                term_description = m.group(1)[:-2] if m.group(1).endswith("tp") else m.group(1)
                files[f"pYPKa_E_{term_description}.gb"] = term.format("gb")
                print(f"pYPKa_E_{term_description}.gb added to file list.")
            else:
                print(f"{term.description} is a linear DNA fragment.")
                term_description =  term.id[:-2] if term.id.endswith("tp") else term.id
                if "pYPKa_ZE_{}.md".format(term_description) not in files:
                    files[term_description+".gb"] = term.format("gb")
                    nbtemp = read_data_file("nb_template_pYPKa_ZE_insert.md")
                    files["pYPKa_ZE_{}.ipynb".format(term_description)] = obj.to_notebook( nbtemp.format(tp=term_description) )
                    pYPKa_clones+="[pYPKa_ZE_{n}](pYPKa_ZE_{n}.ipynb)  \n".format(n=term_description)
                print(f"{tp_vector} added to file list.")

            x = "pYPK0_{}_{}_{}".format(prom_description, gene_description, term_description)

            if pYPKa_A or nbflag:
                nbtemp = read_data_file("nb_template_pYPK0_tp_gene_tp.md")
                files[x+".ipynb"] = obj.to_notebook(nbtemp.format(tpz=prom_description,
                                                    gene=gene_description,
                                                    tpe=term_description))
                print(f"{x}.ipynb added to file list. (pYPKa_A)")
            else:
                nbtemp = read_data_file("nb_template_pYPK0_tp_gene_tp_gap_repair.md")
                files[x+".ipynb"] = obj.to_notebook(nbtemp.format(tpz=prom_description,
                                                    gene=gene.locus,
                                                    tpe=term_description))
                print(f"{x}.ipynb added to file list (gap repair).")

            nbflag=False

            cas_vectors.append(f"{x}.gb")
            tp_gene_tp_links+="[{}]({}.ipynb)  \n".format(x, x)

            promoters.append(prom_description)
            genes.append(gene_description)
            terminators.append(term_description)

    nbtemp = read_data_file("nb_template_pYPK0_pw_from_sequence.md")

    cas_vectors_code = "".join(f'"{v}",' for v in cas_vectors)

    pwname = "pYPK0"

    for x,y,z in zip(promoters,genes,terminators):
        pwname+= f"_{x}_{y}"
    pwname+= f"_{z}"

    print(f"Final pathway name is {pwname}")

    files["pw.ipynb"] = obj.to_notebook( nbtemp.format(name=pwname,
                                         tp_gene_tp_links = tp_gene_tp_links,
                                         cas_vectors=cas_vectors_code,
                                         pYPKa_clones=pYPKa_clones,
                                         length=len(genes)))

    return files, error


def execute_notebooks(files, working_dir="ypkassembly", print=print):

    error = False

    try:
        pathlib.Path(working_dir).mkdir(parents=True)
    except OSError as exception:
        if exception.errno == errno.EEXIST:
            print(f"The {working_dir} directory already exists! Please delete or choose another name.")
        else:
            print(f"The {working_dir} directory could not be created")
        return None

    print(f"created subdirectory {working_dir}")

    for name, content in ((n, c) for n, c in list(files.items()) if not n.endswith(".ipynb")):
        mode = {True:"wb", False:"w"}[hasattr(content, "decode")]
        with open(os.path.join(working_dir, name), mode) as f: f.write(content)
        #del files[name]

    print("saved files and images")

    from nbconvert.preprocessors import CellExecutionError

    ep = ExecutePreprocessor()
    ep.timeout = 1200              # seconds
    ep.interrupt_on_timeout = True

    print("executing notebooks.")

    for name, nb in ((n, c) for n, c in tqdm(files.items()) if n.endswith(".ipynb")):
        newpath = os.path.join(working_dir,os.path.splitext(name)[0]+".ipynb")
        try:
            out, slask = ep.preprocess(nb, {'metadata': {'path': working_dir}})
        except CellExecutionError:
            print(f"Error executing the notebook {name}")
            print(f"See notebook for the traceback.")
        else:
            nbformat.write(out, newpath)
        finally:
            nbformat.write(nb, newpath)

    fl = newpath

    return fl, error

def pathway(text           : Any,
            working_dir    : str ="ypkassembly",
            pYPKa_A        : bool = True,
            print = print):

    files, error = generate(text,
                            pYPKa_A,
                            print=print)

    if error:
        print("errors found. No execution of notebooks.")
    else:
        fl,error = execute_notebooks(files, working_dir, print=print)
    return fl, error



def pYPKa_ZE_from_csv( text           : str,
                       email          : str = "me@example.com",
                       backbone       : str = "pYPKa"):

    error = False
    info = []

    try:
        f = open(text, newline='')
    except IOError:
        pass
    else:
        text = f.read()
    finally:
        f.close()

    dialect = csv.Sniffer().sniff(text)
    reader = csv.reader(text.strip().splitlines(), dialect)

    output = io.StringIO()
    writer = csv.writer(output, quoting=csv.QUOTE_NONNUMERIC)

    for row in reader:
        if len(row)==4:
            info.append(tuple(row))
        elif len(row)==2:
            info.append(tuple(row+[None, None]))
        else:
            error = True
            print("Need two or four items per row.")
            print("geneId, Genbank_accession or")
            print("geneId, Genbank_accession, fwprimer, rvprimer")
            print("|".join(row))
        writer.writerow(row)

    files = { "standard_primers.txt"     : read_data_file("standard_primers.txt"),
              "pYPK_ZE.png"              : read_bin_file("pYPK_ZE.png") }

    files["INDATA.csv"] = output.getvalue()

    files[f"{backbone}.gb"] = read_data_file(f"{backbone}.gb")

    nbtemp_primers = read_data_file(f"nb_template_{backbone}_ZE_genbank.md")
    nbtemp_wo_primers = read_data_file(f"nb_template_{backbone}_ZE_genbank_primer_design.md")

    from pydna.myprimers import primerlist

    list_primers = primerlist()

    obj = notedown.MarkdownReader()

    for insert,gblink,f,r in info:
        print(insert, gblink, f, r)

    for insert,gblink,f,r in info:

        if f and r:
            fp = list_primers[int(f[:3])]
            rp = list_primers[int(r[:3])]

            nb = nbtemp_primers.format( insert=insert,
                                        email=email,
                                        gblink=gblink,
                                        fpn = fp.name,
                                        fps = str(fp.seq),
                                        rpn = rp.name,
                                        rps = str(rp.seq))
        else:
            nb = nbtemp_wo_primers.format( insert=insert,
                                           email=email,
                                           gblink=gblink)

        files[f"{backbone}_ZE_{insert}.ipynb"] = obj.to_notebook(nb)

    return files, error


def pYPKX_A_from_csv( text           : str,
                      email          : str = "me@example.com",
                      backbone       : str = "pYPKa"):

    error = False
    info = []

    try:
        f = open(text, newline='')
    except IOError:
        pass
    else:
        text = f.read()
    finally:
        f.close()

    dialect = csv.Sniffer().sniff(text)
    reader = csv.reader(text.strip().splitlines(), dialect)

    output = io.StringIO()
    writer = csv.writer(output, quoting=csv.QUOTE_NONNUMERIC)

    for row in reader:
        if len(row)==4:
            info.append(tuple(row))
        elif len(row)==2:
            info.append(tuple(row+[None,None]))
        else:
            error = True
            print("Need two or four items per row.")
            print("geneId, Genbank_accession or")
            print("geneId, Genbank_accession, fwprimer, rvprimer")
            print("|".join(row))
        writer.writerow(row)

    files = { "standard_primers.txt"     : read_data_file("standard_primers.txt"),
              "pYPK_A.png"               : read_bin_file("pYPK_A.png") }

    files["INDATA.csv"] = output.getvalue()

    files[f"{backbone}.gb"] = read_data_file(f"{backbone}.gb")

    nbtemp_primers    = read_data_file(f"nb_template_{backbone}_A_genbank.md")
    nbtemp_wo_primers = read_data_file(f"nb_template_{backbone}_A_genbank_primer_design.md")

    from pydna.myprimers import list_primers

    obj = notedown.MarkdownReader()

    for insert, gblink, f, r  in info:
        print(insert, f,r, gblink)

    for insert,gblink,f,r in info:

        if f and r:
            fp = list_primers[int(f[:3])]
            rp = list_primers[int(r[:3])]

            nb = nbtemp_primers.format( insert=insert,
                                        email=email,
                                        gblink=gblink,
                                        fpn = fp.name,
                                        fps = str(fp.seq),
                                        rpn = rp.name,
                                        rps = str(rp.seq))
        else:
            nb = nbtemp_wo_primers.format( insert=insert,
                                           email=email,
                                           gblink=gblink)

        files[f"{backbone}_A_{insert}.ipynb"] = obj.to_notebook(nb)

    return files, error

def main():

    error = False

    appdir = pathlib.Path(appdirs.user_data_dir("ypkpathway")).joinpath("constructs")

    appdir.mkdir(parents=True, exist_ok=True)

    try:
        arguments = docopt.docopt(__doc__)
    except docopt.DocoptExit as e:
        sys.exit(e)


    working_dir = f"ypk_assembly_{datetime.now().isoformat()}"

    if arguments["<dir>"]:
        working_dir = str(arguments["<dir>"])

    if arguments["--version"]:
        from ._version import get_versions
        __version__ = get_versions()["version"][:5]
        del get_versions
        print("ypkpathway version:",__version__)
        print("     pydna version:",pydna.__version__)

    if arguments["assembly_from_file"] and arguments["<path>"]:
        if arguments["--no_pYPKa_A"]:
            pYPKa_A = False
        else:
            pYPKa_A = True
        print("Assembly started from file. (This might take a while...)")
        files, error = generate( str(arguments["<path>"] ),
                                 pYPKa_A=pYPKa_A )

    elif arguments["assembly_from_name"] and arguments["<name>"]:
        print("Assembly started from name. (This might take a while...)")
        files, error = generate(str(arguments["<name>"]))

    elif arguments["pYPKa_A_from_csv"] and arguments["<path>"] and arguments["--email"]:
        print("pYPKa_A clones to ge generated from csv file. (This might take a while...)")
        files, error  = pYPKX_A_from_csv(arguments["<path>"],
                                         email=arguments["--email"],
                                         backbone="pYPKa")

    elif arguments["pYPK0_A_from_csv"] and arguments["<path>"]and arguments["--email"]:
        print("pYPK0_A clones to ge generated from csv file.  (This might take a while...)")
        files, error  = pYPKX_A_from_csv(arguments["<path>"],
                                         email=arguments["--email"],
                                         backbone="pYPK0")

    elif arguments["pYPKa_ZE_from_csv"] and arguments["<path>"]and arguments["--email"]:
        print("pYPKa_ZE clones to ge generated from csv file.  (This might take a while...)")
        files, error  = pYPKa_ZE_from_csv(arguments["<path>"],
                                          email=arguments["--email"],
                                          backbone="pYPKa")

    if not error:

        filepath, error = execute_notebooks(files, working_dir)

        if error:
            print("Notebooks executed with errors")

        print("# If you want to open the notebooks, type:")
        print()
        print(f"    jupyter notebook {filepath}")
        print("\nor\n")
        print(f"    jupyter lab {filepath}  ")


        # #filepath =  os.path.join(working_dir, "pw.ipynb")

        # if platform.system() == 'Darwin':       # macOS
        #     subprocess.run(('open', filepath))
        # elif platform.system() == 'Windows':    # Windows
        #     os.startfile(filepath)
        # else:                                   # linux variants
        #     subprocess.run(('xdg-open', filepath))
    else:
        print("Errors were found in the preparation for assembly")
        print("These should be possible to identify by inspecting")
        print("the log text above.")



if __name__ == "__main__":
    main()


