#!/usr/bin/env python
# -*- coding: utf-8 -*-
u"""Usage: ypkpathway <path>
          ypkpathway -h|--help
          ypkpathway -v|--version

Arguments:
    <path>  path to data file containing sequences to be assembled

Options:
    -h, --help      Show this screen.
    -v, --version   Show version.
"""

import io
import re
#from time import gmtime, strftime
from IPython import nbformat
import notedown
import sys
import os
import errno
import codecs
from IPython.nbconvert.preprocessors.execute import ExecutePreprocessor


from pkg_resources import resource_filename

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

def pathway(pth):

    pYPK0 = pydna.read(read_data_file("pYPK0.gb"))
    pYPKa = pydna.read(read_data_file("pYPKa.gb"))

    from Bio.Restriction import ZraI, AjiI, EcoRV

    files = {"primers.fasta" : read_data_file("primers.fasta"),
             "pYPKa.gb"      : read_data_file("pYPKa.gb"),
             "pYPKpw.gb"     : read_data_file("pYPKpw.gb")}

    cas_vectors = u""
    tp_gene_tp_links = u""

    while pth:
        first = pth.pop(0)
        # is sequence a tp-gene-tp vector?
        if cloned(pYPK0, (ZraI, EcoRV),  first):
            fn = first.description+u".gb"
            files[fn] = first.format("gb")
            cas_vectors+= fn+u"\n"
            tp_gene_tp_links+= u"\n[{}]({})\n".format( first.description, fn )
        else:
            try:
                middle = pth.pop(0)
                last   = pth.pop(0)
            except IndexError:
                raise Exception(u"not enough sequences")

            prom, gene, term = first, middle, last

            if cloned(pYPKa, ZraI,  prom):
                desc_re  = re.compile("pYPKa_Z_([^\d\W]\w{2,15})tp")
            else:
                n = "pYPKa_Z_{}tp".format(prom.id)
                files[prom.id+u".gb"] = prom.format("gb")
                nbtemp = read_data_file("nb_template_pYPKa_ZE_insert.md")
                files[n+u".md"] = nbtemp.format(tp=prom.id)

            if cloned(pYPKa, AjiI,  gene):
                desc_re  = re.compile("pYPKa_A_([^\d\W]\w{2,15})tp")
            else:
                n = u"pYPKa_A_{}".format(gene.id)
                files[gene.id+u".gb"] = gene.format("gb")
                nbtemp = read_data_file("nb_template_pYPKa_A_insert.md")
                files[n+u".md"] = nbtemp.format(insert=gene.id)

            if cloned(pYPKa, EcoRV, term):
                desc_re  = re.compile("pYPKa_E_([^\d\W]\w{2,15})tp")
            else:
                n = u"pYPKa_E_{}tp".format(term.id)
                files[term.id+u".gb"] = term.format("gb")
                nbtemp = read_data_file("nb_template_pYPKa_ZE_insert.md")
                files[n+u".md"] = nbtemp.format(tp=term.id)

            nbtemp = read_data_file("nb_template_pYPK0_tp_gene_tp.md")
            x= "pYPK0_{}tp_{}_{}tp".format(prom.id,gene.id,term.id)
            files[x+u".md"] = nbtemp.format(tpz=prom.id,gene=gene.id,tpe=term.id)

            cas_vectors+=u"\n"+x+u".gb\n"

    nbtemp = read_data_file("nb_template_pYPK0_pw.md")

    pwname = u"pYPK0_"

    pwnb = nbtemp.format(name=pwname,
                         tp_gene_tp_links = tp_gene_tp_links,
                         cas_vectors=add_space(cas_vectors, 17))

    files[u"pw.md"] = pwnb

    return files



def main():
    import docopt
    try:
        arguments = docopt.docopt(__doc__)
    except docopt.DocoptExit as e:
        print e.message

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

        pw_list = pydna.parse( text )

        try:
            os.makedirs("ypk_assembly")
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        pw = pathway( pw_list )

        for name, content in pw.items():
            with codecs.open(os.path.join("ypk_assembly", name), "w", "utf8") as f:
                f.write(content)
        print u"Assembly finished! files written to folder ypkpathway"


if __name__ == "__main__":

    obj = notedown.MarkdownReader()

    with open(u"../tests/pth2.txt", "rU") as f: text = f.read()

    import pydna

    pw = pathway( pydna.parse( text ) )

    print "create subdirectory ypk_assembly"

    try:
        os.makedirs("ypk_assembly")
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise

    os.chdir("ypk_assembly")


    for name, content in ((n, c) for n, c in pw.items() if not n.endswith(".md")):
        print "saving file {}".format(name)
        with open(name,"w") as f: f.write(content)

    for name, content in ((n, c) for n, c in pw.items() if n.endswith(".md")):
        print "saving file {}".format(name)
        newname = os.path.splitext(name)[0]+".ipynb"
        nb = nbformat.writes(obj.to_notebook(content))
        with open(newname,"w") as f: f.write(nb)

    pp = ExecutePreprocessor()
    pp.timeout = 120 # seconds
    pp.interrupt_on_timeout = True

    print
    print "executing pYPKa notebooks"

    for name in (f for f in os.listdir(".") if re.match("pYPKa.+\.ipynb", f)):
        print name
        with io.open(name, 'r', encoding='utf-8') as f: nb = nbformat.read(f, 4)
        nb_executed, resources = pp.preprocess(nb, resources={})
        nbformat.write(nb, name)

    print
    print "executing pYPK0 notebooks"

    for name in (f for f in os.listdir(".") if re.match("pYPK0.+\.ipynb", f)):
        print name
        with io.open(name, 'r', encoding='utf-8') as f: nb = nbformat.read(f, 4)
        nb_executed, resources = pp.preprocess(nb, resources={})
        nbformat.write(nb, name)

    print "executing final pathway notebook"
    with io.open("pw.ipynb", 'r', encoding='utf-8') as f: nb = nbformat.read(f, 4)
    nb_executed, resources = pp.preprocess(nb, resources={})
    nbformat.write(nb, "pw.ipynb")
    print u"Assembly finished! files written to folder ypkpathway"
