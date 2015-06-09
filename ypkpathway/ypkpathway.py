#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""
Usage: ypkpathway <path> [<dir>]
       ypkpathway -h|--help
       ypkpathway -v|--version

Arguments:
    <path>  path to data file containing sequences to be assembled

    <dir>   Directory to put generated sequence files, defaults to
            <ypk_assembly> in the current working directory.

Options:
    -h, --help      Show this screen.
    -v, --version   Show version.
"""

import io
import re
#from time import gmtime, strftime
import sys
import os
import subprocess
import errno
import codecs
import docopt
from IPython import nbformat
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

def pathway(pth, dir_="ypkassembly"):

    pYPK0 = pydna.read(read_data_file("pYPK0.gb"))
    pYPKa = pydna.read(read_data_file("pYPKa.gb"))

    from Bio.Restriction import ZraI, AjiI, EcoRV

    files = {"primers.fasta" : read_data_file("primers.fasta"),
             "pYPKa.gb"      : read_data_file("pYPKa.gb"),
             "pYPKpw.gb"     : read_data_file("pYPKpw.gb")}

    cas_vectors = u""
    tp_gene_tp_links = u""
    pYPKa_clones=u""
    pwname = u"pYPK0"

    while pth:
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
                n = "pYPKa_Z_{}".format(prom.id)
                files[prom.id+u".gb"] = prom.format("gb")
                nbtemp = read_data_file("nb_template_pYPKa_ZE_insert.md")
                files[n+u".md"] = nbtemp.format(tp=prom.id)
                prom_description = prom.id
                pYPKa_clones+=u"[{}]({}.ipynb)  \n".format(n, n)

            if cloned(pYPKa, AjiI,  gene):
                m = re_A.search(gene.description)
                if not m:
                    raise Exception( "{} is a pYPKa_A_gene sequence but was incorrectly named.".format(gene.description))
                gene_description = m.group(1)
                files[m.group(0)+u".gb"] = gene.format("gb")
            else:
                n = u"pYPKa_A_{}".format(gene.locus)
                files[gene.locus+u".gb"] = gene.format("gb")
                nbtemp = read_data_file("nb_template_pYPKa_A_insert.md")
                files[n+u".md"] = nbtemp.format(insert=gene.locus)
                gene_description = gene.locus
                pYPKa_clones+=u"[{}]({}.ipynb)  \n".format(n, n)

            if cloned(pYPKa, EcoRV, term):
                m = re_E.search(term.description)
                if not m:
                    raise Exception( "{} is a pYPKa_A_gene sequence but was incorrectly named.".format(gene.description))
                term_description = m.group(1)
                files[m.group(0)+u".gb"] = term.format("gb")
            else:
                n = u"pYPKa_E_{}".format(term.id)
                files[term.id+u".gb"] = term.format("gb")
                nbtemp = read_data_file("nb_template_pYPKa_ZE_insert.md")
                files[n+u".md"] = nbtemp.format(tp=term.id)
                term_description = term.id
                pYPKa_clones+=u"[{}]({}.ipynb)  \n".format(n, n)

            nbtemp = read_data_file("nb_template_pYPK0_tp_gene_tp.md")
            x = "pYPK0_{}_{}_{}".format(prom_description, gene_description, term_description)
            files[x+u".md"] = nbtemp.format(tpz=prom_description,gene=gene_description,tpe=term_description)
            #pYPK0_clones+=
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

    print "created subdirectory", dir_
    print

    os.chdir(dir_)

    print "saving files.."

    for name, content in ((n, c) for n, c in files.items() if not n.endswith(".md")):
        print name
        with open(name,"w") as f: f.write(content)

    for name, content in ((n, c) for n, c in files.items() if n.endswith(".md")):
        newname = os.path.splitext(name)[0]+".ipynb"
        nb = nbformat.writes(obj.to_notebook(content))
        print newname
        with open(newname,"w") as f: f.write(nb)


    pp = ExecutePreprocessor()
    pp.timeout = 120 # seconds
    pp.interrupt_on_timeout = True


    print
    print "executing pYPKa notebooks.."

    shell = InteractiveShell.instance()
    new_primers = []
    d={}
    for name in (f for f in os.listdir(".") if re.match("pYPKa.+\.ipynb", f)):
        print name
        with io.open(name, 'r', encoding='utf-8') as f: nb = nbformat.read(f, 4)
        nb_executed, resources = pp.preprocess(nb, resources={})
        nbformat.write(nb, name)
        for cell in nb.cells:
            if cell.cell_type == 'code':
                code = shell.input_transformer_manager.transform_cell(cell.source)
                exec( code, d )
        new_primers.extend( (d["fp"], d["rp"]) )
    else:
        print "No pYPKa notebooks found."

    print "executing pYPK0 notebooks.."

    for name in (f for f in os.listdir(".") if re.match("pYPK0.+\.ipynb", f)):
        print name
        with io.open(name, 'r', encoding='utf-8') as f: nb = nbformat.read(f, 4)
        nb_executed, resources = pp.preprocess(nb, resources={})
        nbformat.write(nb, name)
    else:
        print "No pYPK0 notebooks found."

    nbtemp = read_data_file("nb_template_pYPK0_pw.md")

    primer_list = "\n".join( p.format("tab") for p in new_primers )

    pwnb = nbtemp.format(name=pwname,
                         tp_gene_tp_links = tp_gene_tp_links,
                         cas_vectors=add_space(cas_vectors, 17),
                         primer_list=primer_list,
                         pYPKa_clones=pYPKa_clones)

    nb = nbformat.writes(obj.to_notebook(pwnb))
    with open("pw.ipynb", "w") as f: f.write(nb)

    print "executing final pathway notebook"
    print "pw.ipynb"
    with io.open("pw.ipynb", 'r', encoding='utf-8') as f: nb = nbformat.read(f, 4)
    nb_executed, resources = pp.preprocess(nb, resources={})
    nbformat.write(nb, "pw.ipynb")

    os.chdir(cwd)

    fl = FileLink("mypathway/pw.ipynb")

    return fl

def main():

    try:
        arguments = docopt.docopt(__doc__)
    except docopt.DocoptExit as e:
        print e.message
        sys.exit(0)

    dir_ = "ypk_assembly"

    if arguments["<dir>"]:
        dir_= arguments["<dir>"]

    if arguments["--version"]:
        from ._version import get_versions
        __version__ = get_versions()["version"][:5]
        del get_versions
        print u"ypkpathway version:",__version__

    if arguments["<path>"]:
        file_ = arguments["<path>"]
        try:
            with open(file_, "rU") as f: text=f.read()
        except IOError:
            print arguments["<path>"], "could not be opened!"
            sys.exit(1)

        print u"Assembly started! (This might take a while...)"

        fl = pathway( pydna.parse(text), dir_ )

        print u"opening IPython notebook {}".format(fl.path)

        subprocess.Popen(["ipython", "notebook", os.path.join(dir_, "pw.ipynb")])

if __name__ == "__main__":
    pass