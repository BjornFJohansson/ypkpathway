#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
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

from __future__ import print_function
import io
import re
#from time import gmtime, strftime
import sys
import os
#import subprocess
import errno
import codecs
import shutil

import docopt
from IPython import nbformat   #, nbconvert
from IPython.nbconvert.preprocessors.execute import ExecutePreprocessor
from IPython.core.interactiveshell import InteractiveShell
from IPython.display import FileLink

import notedown
import pydna
from pkg_resources import resource_filename

re_cas  = re.compile("pYPK0_([^\d\W]\w{2,15})_([^\d\W]\w{2,15})_([^\d\W]\w{2,15})")
re_cas  = re.compile("pYPK0_([^_]{2,15})_([^_]{2,15})_([^_]{2,15})")
re_Z    = re.compile("pYPKa_Z_([^\d\W]\w{2,15})")
re_A    = re.compile("pYPKa_A_([^\d\W]\w{2,15})")
re_E    = re.compile("pYPKa_E_([^\d\W]\w{2,15})")

def add_space(s, n):
    return  u"\n".join(u"{}{}".format(u" "*n, line) for line in s.splitlines())

def cloned(vector, enzyme, candidate):
    if len(candidate) <= len(vector):
        return 0
    candidate2 = str(candidate.seq.tolinear()*2).lower()
    linear_vector = vector.cut(enzyme).pop(0)
    if str(linear_vector.seq).lower() in candidate2:
        return len(candidate) - len(vector)
    return 0

def read_data_file(name):
    with codecs.open( resource_filename("ypkpathway", os.path.join("data", name)), "r", "utf-8") as f: data = f.read()
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

    log=u""

    pYPK0 = pydna.read(read_data_file("pYPK0.gb"))
    pYPKa = pydna.read(read_data_file("pYPKa.gb"))

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

    cas_vectors = u""
    tp_gene_tp_links = u""
    pYPKa_clones=u""
    pwname = u"pYPK0"
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
            fn = first.description+u".gb"
            files[fn] = first.format("gb")
            cas_vectors+= fn+u"\n"
            tp_gene_tp_links+= u"\n[{}]({})\n".format( first.description, fn )
            tp1_description  = m.group(1)
            gene_description = m.group(2)
            tp2_description  = m.group(3)
            genes+=1
        else:
            try:
                middle = pth.pop(0)
                last   = pth.pop(0)
            except IndexError:
                raise Exception(u"not enough sequences")

            prom, gene, term = first, middle, last

            if cloned(pYPKa, ZraI,  prom):
                m = re_Z.search(prom.description)
                if not m:
                    raise Exception( "{} is a pYPKa_A_gene sequence but was incorrectly named.".format(gene.description))
                prom_description = m.group(1)
                files[m.group(0)+u".gb"] = prom.format("gb")
            else:
                #print("Z"+str(files.has_key("pYPKa_ZE_{}.md".format(prom.id)))+prom.id)
                if not files.has_key(u"pYPKa_ZE_{}.md".format(prom.id)):
                    files[prom.id+u".gb"] = prom.format("gb")
                    nbtemp = read_data_file("nb_template_pYPKa_ZE_insert.md")
                    files[u"pYPKa_ZE_{}.md".format(prom.id)] = nbtemp.format(tp=prom.id)
                    pYPKa_clones+=u"[pYPKa_ZE_{n}](pYPKa_ZE_{n}.ipynb)  \n".format(n=prom.id)
                prom_description = prom.id

            if cloned(pYPKa, AjiI,  gene):
                m = re_A.search(gene.description)
                if not m:
                    raise Exception( "{} is a pYPKa_A_gene sequence but was incorrectly named.".format(gene.description))
                gene_description = m.group(1)
                files[m.group(0)+u".gb"] = gene.format("gb")
                if not pYPKa_A:
                    nbflag=True

            else:
                n = u"pYPKa_A_{}".format(gene.locus)
                files[gene.locus+u".gb"] = gene.format("gb")
                if pYPKa_A:
                    nbtemp = read_data_file("nb_template_pYPKa_A_insert.md")
                    files[n+u".md"] = nbtemp.format(insert=gene.locus)
                    gene_description = gene.locus
                    pYPKa_clones+=u"[{}]({}.ipynb)  \n".format(n, n)
                else:
                    gene_description = gene.locus

            if cloned(pYPKa, EcoRV, term):
                m = re_E.search(term.description)
                if not m:
                    raise Exception( "{} is a pYPKa_A_gene sequence but was incorrectly named.".format(gene.description))
                term_description = m.group(1)
                files[m.group(0)+u".gb"] = term.format("gb")
            else:
                #print("E"+str(files.has_key("pYPKa_ZE_{}.md".format(term.id)))+term.id)
                if not files.has_key(u"pYPKa_ZE_{}.md".format(term.id)):
                    files[term.id+u".gb"] = term.format("gb")
                    nbtemp = read_data_file("nb_template_pYPKa_ZE_insert.md")
                    files[u"pYPKa_ZE_{}.md".format(term.id)] = nbtemp.format(tp=term.id)
                    pYPKa_clones+=u"[pYPKa_ZE_{n}](pYPKa_ZE_{n}.ipynb)  \n".format(n=term.id)
                term_description = term.id

            x = "pYPK0_{}_{}_{}".format(prom_description, gene_description, term_description)

            if pYPKa_A or nbflag:
                nbtemp = read_data_file("nb_template_pYPK0_tp_gene_tp.md")
                files[x+u".md"] = nbtemp.format(tpz=prom_description,
                                                gene=gene_description,
                                                tpe=term_description)
            else:
                nbtemp = read_data_file("nb_template_pYPK0_tp_gene_tp_gap_repair.md")
                files[x+u".md"] = nbtemp.format(tpz=prom_description,
                                                gene=gene.locus,
                                                tpe=term_description)
            nbflag=False

            cas_vectors+=u"\n"+x+u".gb\n"
            tp_gene_tp_links+=u"[{}]({}.ipynb)  \n".format(x, x)




        pwname+="_{}".format(gene_description)

    ###########################################################################

    obj = notedown.MarkdownReader()

    cwd = os.getcwd()

    try:
        os.makedirs(dir_)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    msg = u"created subdirectory {}\n".format(dir_)
    print(msg)
    log+=msg

    os.chdir(dir_)

    msg = u"\nsaving files sequence files and images..\n"
    print(msg)
    log+=msg

    for name, content in sorted((n, c) for n, c in files.items() if not n.endswith(".md")):
        msg = u"\nsaving: "+name
        print(msg)
        log+=msg
        with open(name,"wb") as f: f.write(content)

    print("\n")
    log+="\n"

    msg = u"\nsaving notebook files ..\n"
    print(msg)
    log+=msg

    for name, content in sorted((n, c) for n, c in files.items() if n.endswith(".md")):
        newname = os.path.splitext(name)[0]+".ipynb"
        msg = u"\nsaving: "+newname
        print(msg)
        log+=msg
        nb = nbformat.write(obj.to_notebook(content), newname)

    pp = ExecutePreprocessor()
    pp.timeout = 120 # seconds
    pp.interrupt_on_timeout = True

    print("\n")
    log+="\n"

    msg = u"\nexecuting pYPKa notebooks..\n"
    print(msg)
    log+=msg

    shell = InteractiveShell.instance()
    #new_primers = []

    g={}
    l={}

    pypkanbs = sorted([f for f in os.listdir(".") if re.match("pYPKa.+\.ipynb", f)])

    if pypkanbs:
        for name in pypkanbs:
            msg = u"\nexecuting: "+name
            print(msg)
            log+=msg
            with io.open(name, 'r', encoding='utf-8') as f: nb = nbformat.read(f, 4)
            nb_executed, resources = pp.preprocess(nb, resources={})
            for cell in nb.cells:
                if cell.cell_type == 'code':
                    code = shell.input_transformer_manager.transform_cell(cell.source)
                    exec code in g, l
            #new_primers.extend( (l["fp"], l["rp"]) )
            nbformat.write(nb, name)
            g={}
            l={}
    else:
        msg = u"\nNo pYPKa notebooks found.\n"
        print(msg)
        log+=msg
    print("\n")
    log+="\n"
    msg = u"\nexecuting pYPK0 notebooks..\n"
    print(msg)
    log+=msg

    g={}
    l={}
    resources={}

    pypk0nbs = sorted([f for f in os.listdir(".") if re.match("pYPK0.+\.ipynb", f)])

    if pypk0nbs:
        for name in pypk0nbs:
            msg = u"\nexecuting: "+name
            print(msg)
            log+=msg
            with io.open(name, 'r', encoding='utf-8') as f: nb = nbformat.read(f, 4)
            nb_executed, resources = pp.preprocess(nb, resources={})
            nbformat.write(nb, name)
            for cell in nb.cells:
                if cell.cell_type == 'code':
                    code = shell.input_transformer_manager.transform_cell(cell.source)
                    exec code in g, l
            #try:
                #new_primers.extend( (l["fp"], l["rp"]) )
            #except KeyError:
            #    pass
            g={}
            l={}
    else:
        msg = u"\nNo pYPK0 notebooks found.\n"
        print(msg)
        log+=msg
    nbtemp = read_data_file("nb_template_pYPK0_pw.md")

    #primer_list = "\n".join( p.format("tab") for p in new_primers )

    #if new_primers:
    #    msg = u"\n\nsaving new_primers.txt..\n"
    #with open("new_primers.txt","wb") as f: f.write("\n".join( p.format("fasta") for p in new_primers ))

    pwnb = nbtemp.format(name=pwname,
                         filename=os.path.basename(dir_),
                         tp_gene_tp_links = tp_gene_tp_links,
                         cas_vectors=add_space(cas_vectors, 17),
                         pYPKa_clones=pYPKa_clones,
                         length=genes)

    nb = nbformat.write(obj.to_notebook(pwnb), "pw.ipynb")

    #nb = nbformat.writes("pw.ipynb", obj.to_notebook(pwnb))
    #with open("pw.ipynb", "w") as f: f.write(nb)

    msg = u"\n\nexecuting final pathway notebook..\n"
    print(msg)
    log+=msg
    msg = u"\nexecuting: pw.ipynb"
    print(msg)
    log+=msg
    with io.open("pw.ipynb", 'r', encoding='utf-8') as f: nb = nbformat.read(f, 4)
    nb_executed, resources = pp.preprocess(nb, resources={})
    nbformat.write(nb, "pw.ipynb")

    #for nb_ in [f for f in os.listdir(".") if f.endswith(".ipynb")]:
    #    subprocess.Popen(["ipython", "nbconvert", os.path.join(dir_, nb_)])

    os.chdir(cwd)

    fl = FileLink(os.path.join(dir_, u"pw.ipynb"))

    #   pp = None

    return fl, log

def main():

    try:
        arguments = docopt.docopt(__doc__)
    except docopt.DocoptExit as e:
        print(e.message)
        sys.exit(0)

    dir_ = u"ypk_assembly"

    if arguments["<dir>"]:
        dir_= unicode(arguments["<dir>"])

    if arguments["--no_pYPKa_A"]:
        pYPKa_A = False
    else:
        pYPKa_A = True

    if arguments["--version"]:
        from _version import get_versions
        __version__ = get_versions()["version"][:5]
        del get_versions
        print(u"ypkpathway version:",__version__)
        print(u"     pydna version:",pydna.__version__)

    if arguments["<path>"]:
        file_ = unicode(arguments["<path>"])
        try:
            with codecs.open(file_, "rU", 'utf8') as f: text=f.read()
        except IOError:
            print(arguments["<path>"], "could not be opened!")
            sys.exit(1)

        #dir_ = os.path.splitext(os.path.basename(file_))[0]
        dir_, ext = os.path.splitext(os.path.abspath(file_))

        print(u"Assembly started! (This might take a while...)")
        #print(file_)
        #print(dir_)
        #print(os.path.abspath(file_))
        #print(os.path.abspath(dir_))
        #print(os.path.splitext(os.path.abspath(file_)))
        #print(os.path.basename(dir_))
        #import sys;sys.exit(42)

        fl, log = pathway( pydna.parse(text), dir_, pYPKa_A=pYPKa_A )

        with codecs.open(os.path.join(dir_, u"log.txt"),"wU","utf8") as f: f.write(log)

        filename = os.path.basename(file_)

        shutil.copy2( filename, os.path.join(dir_, u"INDATA_"+os.path.basename((dir_)+u".txt")))

        print(u"opening IPython notebook {}".format(fl.path))

        #subprocess.Popen(["ipython", "notebook", os.path.join(dir_, "pw.ipynb")])

def pathway_(x,y, print=print):
    print("abc")
    fl = FileLink(os.path.join("dir", "pw.ipynb"))
    return fl


if __name__ == "__main__":
    main()
