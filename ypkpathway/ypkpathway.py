#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
ypkpathway_cli

Usage: ypkpathway_cli <path> [<dir>] [--no_pYPKa_A]
       ypkpathway_cli -h|--help
       ypkpathway_cli -v|--version

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

from pkg_resources import resource_filename

import pydna
from pydna.readers import read
from pydna.parsers import parse


import nbformat
from nbconvert.preprocessors.execute import ExecutePreprocessor

from IPython.core.interactiveshell import InteractiveShell
from IPython.display import FileLink, FileLinks

import notedown

re_cas  = re.compile("pYPK0_([^\d\W]\w{2,15})_([^\d\W]\w{2,15})_([^\d\W]\w{2,15})")
re_cas  = re.compile("pYPK0_([^_]{2,15})_([^_]{2,15})_([^_]{2,15})")
re_Z    = re.compile("pYPKa_Z_([^\d\W]\w{2,15})")
re_A    = re.compile("pYPKa_A_([^\d\W]\w{2,15})")
re_E    = re.compile("pYPKa_E_([^\d\W]\w{2,15})")

def add_space(s, n):
    return  "\n".join("{}{}".format(" "*n, line) for line in s.splitlines())

def cloned(vector, enzyme, candidate):
    if len(candidate) <= len(vector):
        return 0
    candidate2 = str(candidate.seq.tolinear()*2).lower()
    linear_vector = vector.cut(enzyme).pop(0)
    if str(linear_vector.seq).lower() in candidate2:
        return len(candidate) - len(vector)
    return 0

def read_data_file(name):
    with open( resource_filename("ypkpathway", os.path.join("data", name)), "r") as f: data = f.read()
    return data

def read_bin_file(name):
    with open( resource_filename("ypkpathway", os.path.join("data", name)), "rb") as f: data = f.read()
    return data

def pathway(pth, dir_="ypkassembly", pYPKa_A=True, print=print):

    if len(pth)==0: # pth has to contain some sequences
        print("No of sequences found.")
        return None, None

    names = [s.name for s in pth] # sequence names has to be unique

    #if len(names)>len(set(names)):
    #    print("Gene names are not unique. Please rename sequences so that each sequence has a unique name.\n")
    #    print("Gene names parsed from Data page:\n\n")
    #    for name in names:
    #        print(name)
    #    return None, None

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
                raise Exception( "{} is a pYPK0 tp-gene_tp sequence but was not correctly named.".format(last.description))
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

    ###########################################################################

    obj = notedown.MarkdownReader()

    cwd = os.getcwd()

    try:
        os.makedirs(dir_)
    except OSError as exception:
        if exception.errno == errno.EEXIST:
            print("The {} directory already exists! Please delete or choose another name.".format(dir_))
        else:
            print("The {} directory could not be created".format(dir_))
        return None, None

    msg = "created subdirectory {}\n".format(dir_)
    print(msg)
    log+=msg

    os.chdir(dir_)

    msg = "\nsaving files sequence files and images..\n"
    print(msg)
    log+=msg

    for name, content in sorted((n, c) for n, c in list(files.items()) if not n.endswith(".md")):
        msg = "\nsaving: "+name
        print(msg)
        log+=msg
        mode = {True:"wb", False:"w"}[hasattr(content, "decode")]
        with open(name, mode) as f:  #with open(name,"wb") as f: 
            f.write(content) 

    print("\n")
    log+="\n"

    msg = "\nsaving notebook files ..\n"
    print(msg)
    log+=msg

    for name, content in sorted((n, c) for n, c in list(files.items()) if n.endswith(".md")):
        newname = os.path.splitext(name)[0]+".ipynb"
        msg = "\nsaving: "+newname
        print(msg)
        log+=msg
        nb = nbformat.write(obj.to_notebook(content), newname)

    pp = ExecutePreprocessor()
    pp.timeout = 120 # seconds
    pp.interrupt_on_timeout = True

    print("\n")
    log+="\n"

    msg = "\nexecuting pYPKa notebooks..\n"
    print(msg)
    log+=msg

    shell = InteractiveShell.instance()
    #new_primers = []

    g={}
    l={}

    pypkanbs = sorted([f for f in os.listdir(".") if re.match("pYPKa.+\.ipynb", f)])

    if pypkanbs:
        for name in pypkanbs:
            msg = "\nexecuting: "+name
            print(msg)
            log+=msg
            with io.open(name, 'r', encoding='utf-8') as f: nb = nbformat.read(f, 4)
            nb_executed, resources = pp.preprocess(nb, resources={})
            for cell in nb.cells:
                if cell.cell_type == 'code':
                    code = shell.input_transformer_manager.transform_cell(cell.source)
                    exec(code, g, l)
            #new_primers.extend( (l["fp"], l["rp"]) )
            nbformat.write(nb, name)
            g={}
            l={}
    else:
        msg = "\nNo pYPKa notebooks found.\n"
        print(msg)
        log+=msg
    print("\n")
    log+="\n"
    msg = "\nexecuting pYPK0 notebooks..\n"
    print(msg)
    log+=msg

    g={}
    l={}
    resources={}

    pypk0nbs = sorted([f for f in os.listdir(".") if re.match("pYPK0.+\.ipynb", f)])

    if pypk0nbs:
        for name in pypk0nbs:
            msg = "\nexecuting: "+name
            print(msg)
            log+=msg
            with io.open(name, 'r', encoding='utf-8') as f: nb = nbformat.read(f, 4)
            nb_executed, resources = pp.preprocess(nb, resources={})
            nbformat.write(nb, name)
            for cell in nb.cells:
                if cell.cell_type == 'code':
                    code = shell.input_transformer_manager.transform_cell(cell.source)
                    exec(code, g, l)
            #try:
                #new_primers.extend( (l["fp"], l["rp"]) )
            #except KeyError:
            #    pass
            g={}
            l={}
    else:
        msg = "\nNo pYPK0 notebooks found.\n"
        print(msg)
        log+=msg
    nbtemp = read_data_file("nb_template_pYPK0_pw.md")

    #primer_list = "\n".join( p.format("tab") for p in new_primers )

    #if new_primers:
    #    msg = u"\n\nsaving new_primers.txt..\n"
    #with open("new_primers.txt","wb") as f: f.write("\n".join( p.format("fasta") for p in new_primers ))

    #print("qwerty")
    #print(pwname)
    #print(os.path.basename(dir_))
    #print(tp_gene_tp_links)
    #print(add_space(cas_vectors, 17))
    #print(pYPKa_clones)
    #print(str(genes))
    #print("123456789")

    pwnb = nbtemp.format(name=pwname,
                         filename=os.path.basename(dir_),
                         tp_gene_tp_links = tp_gene_tp_links,
                         cas_vectors=add_space(cas_vectors, 17),
                         pYPKa_clones=pYPKa_clones,
                         length=genes)

    nb = nbformat.write(obj.to_notebook(pwnb), "pw.ipynb")

    #nb = nbformat.writes("pw.ipynb", obj.to_notebook(pwnb))
    #with open("pw.ipynb", "w") as f: f.write(nb)

    msg = "\n\nexecuting final pathway notebook..\n"
    print(msg)
    log+=msg
    msg = "\nexecuting: pw.ipynb"
    print(msg)
    log+=msg
    with io.open("pw.ipynb", 'r', encoding='utf-8') as f: nb = nbformat.read(f, 4)
    nb_executed, resources = pp.preprocess(nb, resources={})
    nbformat.write(nb, "pw.ipynb")

    #for nb_ in [f for f in os.listdir(".") if f.endswith(".ipynb")]:
    #    subprocess.Popen(["ipython", "nbconvert", os.path.join(dir_, nb_)])

    os.chdir(cwd)

    fl = FileLink(os.path.join(dir_, "pw.ipynb"))

    #   pp = None

    return fl, log

def main():

    try:
        arguments = docopt.docopt(__doc__)
    except docopt.DocoptExit as e:
        #print(e)
        sys.exit(0)

    dir_ = "ypk_assembly"

    if arguments["<dir>"]:
        dir_= str(arguments["<dir>"])

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
        file_ = str(arguments["<path>"])
        try:
            with open(file_, "r") as f: text=f.read()
        except IOError:
            print(arguments["<path>"], "could not be opened!")
            sys.exit(1)

        #dir_ = os.path.splitext(os.path.basename(file_))[0]
        dir_, ext = os.path.splitext(os.path.abspath(file_))

        print("Assembly started! (This might take a while...)")
        #print(file_)
        #print(dir_)
        #print(os.path.abspath(file_))
        #print(os.path.abspath(dir_))
        #print(os.path.splitext(os.path.abspath(file_)))
        #print(os.path.basename(dir_))
        #import sys;sys.exit(42)

        fl, log = pathway( parse(text), dir_, pYPKa_A=pYPKa_A )

        with open(os.path.join(dir_, "log.txt"),"w") as f: f.write(log)

        filename = os.path.basename(file_)

        shutil.copy2( filename, os.path.join(dir_, "INDATA_"+os.path.basename((dir_)+".txt")))

        print("opening IPython notebook {}".format(fl.path))

        subprocess.Popen(["jupyter", "notebook", os.path.join(dir_, "pw.ipynb")])

def pathway_(x,y, print=print):
    print("abc")
    fl = FileLink(os.path.join("dir", "pw.ipynb"))
    return fl

def pYPKa_ZE_ipynb_generator(tp, dir_="pYPKa_ZE_vectors"):

    cwd = os.getcwd()

    try:
        os.makedirs(dir_)
    except OSError as exception:
        if exception.errno == errno.EEXIST:
            pass
        else:
            print("The {} directory could not be created".format(dir_))
            return None

    os.chdir(dir_)

    with open("standard_primers.txt","w") as f: f.write(read_data_file("standard_primers.txt"))
    with open("pYPKa.gb","w") as f: f.write(read_data_file("pYPKa.gb"))
    with open("pYPK_ZE.png","w") as f: f.write(read_bin_file("pYPK_ZE.png"))
    with open(tp.id+".gb","w") as f: f.write(tp.format("gb"))

    nbtemp = read_data_file("nb_template_pYPKa_ZE_insert.md")

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

    from io import StringIO
    old_stdout = sys.stdout
    redirected_output = sys.stdout = StringIO()

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


