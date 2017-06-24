#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Usage: ypkpathway <path>
          ypkpathway -h|--help
          ypkpathway -v|--version

Arguments:
    <path>  path to data file containing sequences to be assembled

Options:
    -h, --help      Show this screen.
    -v, --version   Show version.
"""



import re
from time import gmtime, strftime
import zipfile
import cStringIO
import sys
import os
import errno
import codecs

from docutils.core 				import publish_string
from docutils.writers.html4css1    import Writer as HisWriter
from pkg_resources 				import resource_filename
from Bio.Restriction 			import ZraI, AjiI, EcoRV

import pydna

pYPKa  = pydna.read( resource_filename('ypkpathway', os.path.join('data', 'pYPKa.txt')))
pYPK0  = pydna.read( resource_filename('ypkpathway', os.path.join('data', 'pYPK0.txt')))
pYPKpw = pydna.read( resource_filename('ypkpathway', os.path.join('data', 'pYPKpw.txt')))


(p577,
 p578,
 p468,
 p467,
 p567,
 p568,
 p775,
 p778,
 p342) = pydna.parse( u'''  >577
                            gttctgatcctcgagcatcttaagaattc
                            >578
                            gttcttgtctcattgccacattcataagt
                            >468
                            gtcgaggaacgccaggttgcccact
                            >467
                            ATTTAAatcctgatgcgtttgtctgcacaga
                            >567
                            GTcggctgcaggtcactagtgag
                            >568
                            GTGCcatctgtgcagacaaacg
                            >775
                            gcggccgctgacTTAAAT
                            >778
                            ggtaaatccggatTAATTAA
                            >342
                            CCTTTTTACGGTTCCTGGCCT''', ds=False)




def add_line_block(s):
    return  u"\n".join( u"|{}".format(line) for line in s.splitlines())

def add_space(s):
    return  u"\n".join( u" {}".format(line) for line in s.splitlines())

def cloned(vector, enzyme, candidate):
    if len(candidate) <= len(vector):
        return 0
    candidate2 = str(candidate.seq.tolinear()*2).lower()
    linear_vector = vector.cut(enzyme).pop(0)
    if str(linear_vector.seq).lower() in candidate2:
        return len(candidate) - len(vector)
    return 0

class pYPKa_clone(object):

    with codecs.open( resource_filename('ypkpathway',os.path.join('data', 'template_pYPKa_X_insert.txt')), "r", "utf-8" ) as f:
        codetempl = f.read()

    with codecs.open( resource_filename('ypkpathway',os.path.join('data', 'template_pYPKa_plan.txt')), "r", "utf-8" ) as f:
        plantempl = f.read()

    minlength   = 16
    maxlength   = 29
    target_tm   = 60

    def __init__(self, enzyme, data):
        fp_tail = "ttaaat"
        rp_tail = "taattaa"
        if enzyme == ZraI:
            fp = p577
            rp = p342
            desc_re  = re.compile("pYPKa_Z_([^\d\W]\w{2,15})tp")
            self.letter = "Z"
        elif enzyme == AjiI:
            fp = p468
            rp = p342
            fp_tail = "aa"
            rp_tail = ""
            desc_re  = re.compile("pYPKa_A_([^\d\W]\w{2,15})")
            self.letter = "A"
        elif enzyme == EcoRV:
            fp = p568
            rp = p342
            desc_re  = re.compile("pYPKa_E_([^\d\W]\w{2,15})tp")
            self.letter = "E"
        else:
            raise Exception(u"Enzyme has to be ZraI, Ajii or EcoRV, but got {}".format(enzyme))

        self.enzyme  = enzyme

        self.insert_length = cloned(pYPKa, enzyme, data)

        if self.insert_length:
            m = desc_re.search(data.description)
            if not m:
                raise Exception(u"{} is a pYPKa_{} sequence but was not correctly named.".format(data.description,
                                                                                                self.letter))
            self.insert_description = m.group(1)
            self.name = m.group(0)
            self.code = "{} = read('{}.txt')".format({ "Z":"first",
                                                       "A":"middle",
                                                       "E":"last"}[self.letter],
                                                      self.name)
            self.files = { "{}.txt".format(self.name)       : data.format("gb").encode("utf-8"),
                           "{}_plan.rst".format(self.name) :  u"This vector was given!"}
            self.rec = data
            self.flag = True

        elif data.linear:
            self.insert_description = data.name # .id .description
            f,r = pydna.cloning_primers(data,
                                        minlength=pYPKa_clone.minlength,
                                        maxlength=pYPKa_clone.maxlength,
                                        fp_tail=fp_tail,
                                        rp_tail=rp_tail)


            self.insert = pydna.pcr(f, r, data)

            self.insert_length = len(self.insert)

            self.insert.description+= "_prd"

            self.name = u"pYPKa_{}_{}{}".format(self.letter,
                                                data.name.split("tp")[0],
                                                {ZraI : "tp",
                                                 AjiI : "",
                                                 EcoRV: "tp" }[enzyme])

            code = pYPKa_clone.codetempl.format( enz = self.enzyme,
                                                 tp  = data.name,
                                                 f   = f.format("fasta"),
                                                 r   = r.format("fasta"),
                                                 vn  = self.name)



            self.rec =    (pYPKa.linearize(self.enzyme) + self.insert).looped().synced("tcgcgcgtttcggtgatgacggtgaaaacctctg")
            self.rec_rv = (pYPKa.linearize(self.enzyme) + self.insert.rc()).looped().synced("tcgcgcgtttcggtgatgacggtgaaaacctctg")

            self.rec.id = self.name[:15]

            self.rec.description = self.name

            plan = pYPKa_clone.plantempl.format( template =  u"{}_template".format(data.name),
                                                 rec = self.name,
                                                 fwd = self.insert.forward_primer.name,
                                                 rev = self.insert.reverse_primer.name,
                                                 figure = add_space(self.insert.figure()),
                                                 program = add_space(self.insert.program()),
                                                 pcr_product =self.insert_description,
                                                 length = len(self.insert),
                                                 enz=enzyme,
                                                 name = self.name,
                                                 line =  "="*len(self.name),
                                                 fp = fp.name,
                                                 rp = rp.name,
                                                 f  = f.name,
                                                 correct_products =  ", ".join(str(f) for f in [len(p) for p in pydna.Anneal((fp,rp,f), self.rec).products]),
                                                 reversed_products =  ", ".join(str(f) for f in [len(p) for p in pydna.Anneal((fp,rp,f), self.rec_rv).products]),
                                                 clone_empty =  ", ".join(str(f) for f in [len(p) for p in pydna.Anneal((fp,rp,f), pYPKa).products]))

            self.files = { "{}_template.txt".format(data.name) : data.format("gb").decode("utf-8"),
                           "{}.txt".format(self.insert_description) : self.insert.format("gb").decode("utf-8"),
                           "{}.txt".format(self.name)               : self.rec.format("gb").decode("utf-8"),
                           "{}.py".format(self.name)                : code,
                           "{}_plan.rst".format(self.name)          : plan}

            self.code = "from {0} import {0} as {1}".format(self.name,
                                                            {ZraI : "first",
                                                             AjiI : "middle",
                                                             EcoRV: "last" }[enzyme])

            self.flag = False
        else:
            raise Exception("{} has to be a linear DNA fragment or a pYPKa_{} clone!".format(data.description,
                                                                                             self.letter))

    def __repr__(self):
        return  "{} {}".format(self.name, {True :  "given", False:  "cloned" }[self.flag])










class pYPK0_tp_gene_tp(object):

    desc_re  = re.compile("pYPK0_([^\d\W]\w{2,15})tp_([^\d\W]\w{2,15})_([^\d\W]\w{2,15})tp")

    with codecs.open(resource_filename('ypkpathway',os.path.join('data','template_pYPK0_tp_gene_tp.txt')), "r", "utf-8" ) as f:
        codetempl = f.read()

    with codecs.open(resource_filename('ypkpathway',os.path.join('data','template_pYPK0_plan.txt')), "r", "utf-8" ) as f:
        plantempl = f.read()

    #pYPK0_E_Z, stuffer = pYPK0.cut((EcoRV, ZraI))
    pYPKpw_lin = pYPKpw.linearize(EcoRV)

    def __init__(self, *args):

        self.pYPKa_ZraI_tp1  = None
        self.pYPKa_AjiI_gene = None
        self.pYPKa_EcoRV_tp2 = None

        if len(args)==1:
            self.seq = args[0]
            m = pYPK0_tp_gene_tp.desc_re.search(self.seq.description)
            if not m:
                raise Exception( "{} is a pYPK0 tp-gene_tp sequence but was not correctly named.".format(last.description))

            self.tp1_description  = m.group(1)
            self.gene_description = m.group(2)
            self.tp2_description  = m.group(3)
            self.pYPKa_clones = []

            self.seq.description = self.seq.description.replace(" ","_")

            self.files = {  "{}.txt".format(self.seq.description) : self.seq.format("gb"),
                            "{}_plan.rst".format(self.seq.description) :  u"This vector was given!"}
            self.code =  "{i} = read('{f}')".format(f= "{}.txt".format(self.seq.description), i= "{}")

        elif len(args)==3:


            self.pYPKa_clones = [pYPKa_clone(ZraI, args[0]),
                                 pYPKa_clone(AjiI, args[1]),
                                 pYPKa_clone(EcoRV,args[2])]

            first  = pydna.pcr( p577, p567, self.pYPKa_clones[0].rec)
            middle = pydna.pcr( p468, p467, self.pYPKa_clones[1].rec)
            last   = pydna.pcr( p568, p578, self.pYPKa_clones[2].rec)

            self.assembly = pydna.Assembly([pYPK0_tp_gene_tp.pYPKpw_lin, first, middle, last], limit=31)
            self.seq = self.assembly.circular_products[0]

            self.tp1_description  = self.pYPKa_clones[0].insert_description
            self.gene_description = self.pYPKa_clones[1].insert_description
            self.tp2_description  = self.pYPKa_clones[2].insert_description

            self.seq.name = "pYPK0_tp_gene_tp"
            self.seq.id = "-"
            #self.seq.annotations['source']   = "123"       # SOURCE
            #self.seq.annotations['organism'] = "123"       # ORGANISM
            #self.seq.annotations['comment']  = "123"       # COMMENT

            self.seq.description =  "pYPK0_{}tp_{}_{}tp".format(self.pYPKa_clones[0].insert_description.split("tp")[0],
                                                                self.pYPKa_clones[1].insert_description,
                                                                self.pYPKa_clones[2].insert_description.split("tp")[0])



            tp_gene_size = [len(p) for p in pydna.Anneal((p577,p467), self.seq).products]
            gene_tp_size = [len(p) for p in pydna.Anneal((p468,p578), self.seq).products]

            plan = pYPK0_tp_gene_tp.plantempl.format( name = self.seq.description,
                                                      tp1  = first.name,
                                                      gene = middle.name,
                                                      tp2  = last.name,
                                                      figure  = add_space(self.seq.small_fig()),
                                                      pcr1=add_space(first.figure()),
                                                      pcr2=add_space(middle.figure()),
                                                      pcr3=add_space(last.figure()),
                                                      prg1=add_space(first.program()),
                                                      prg2=add_space(middle.program()),
                                                      prg3=add_space(last.program()),
                                                      tp1_name  = self.pYPKa_clones[0].name+"_pcr_prd",
                                                      gene_name = self.pYPKa_clones[1].name+"_pcr_prd",
                                                      tp2_name  = self.pYPKa_clones[2].name+"_pcr_prd",
                                                      tmp1 = self.pYPKa_clones[0].name,
                                                      tmp2 = self.pYPKa_clones[1].name,
                                                      tmp3 = self.pYPKa_clones[2].name,
                                                      line =  "="*len(self.seq.description),
                                                      p1 = first.forward_primer.name,
                                                      p2 = first.reverse_primer.name,
                                                      p3 = middle.forward_primer.name,
                                                      p4 = middle.reverse_primer.name,
                                                      p5 = last.forward_primer.name,
                                                      p6 = last.reverse_primer.name,
                                                      correct_first_tp_gene_prd =  ", ".join(str(s) for s in tp_gene_size),
                                                      missing_first_tp_prd =  ", ".join(str(s-self.pYPKa_clones[0].insert_length) for s in tp_gene_size),
                                                      missing_gene_prd1 =  ", ".join(str(s-self.pYPKa_clones[1].insert_length) for s in tp_gene_size),
                                                      empty_prd1 =  ", ".join(str(s-self.pYPKa_clones[0].insert_length-self.pYPKa_clones[1].insert_length   ) for s in tp_gene_size),
                                                      correct_gene_tp_prd =  ", ".join(str(s) for s in gene_tp_size),
                                                      missing_gene_prd2 =  ", ".join(str(s-self.pYPKa_clones[1].insert_length) for s in gene_tp_size),
                                                      missing_last_tp_prd =  ", ".join(str(s-self.pYPKa_clones[2].insert_length) for s in gene_tp_size),
                                                      empty_prd2 =  ", ".join(str(s-self.pYPKa_clones[1].insert_length-self.pYPKa_clones[2].insert_length) for s in gene_tp_size),)


            code = [ self.pYPKa_clones[0].code,
                     self.pYPKa_clones[1].code,
                     self.pYPKa_clones[2].code]

            code = pYPK0_tp_gene_tp.codetempl.format(code =  "\n".join(code),
                                                     vn   =  self.seq.description)

            self.files = {self.pYPKa_clones[0].name+"_pcr_prd.txt"  : first.format("gb").decode("utf-8"),
                          self.pYPKa_clones[1].name+"_pcr_prd.txt"  : middle.format("gb").decode("utf-8"),
                          self.pYPKa_clones[2].name+"_pcr_prd.txt"  : last.format("gb").decode("utf-8"),


                           "{}.py".format(self.seq.description)  : code,
                           "{}.txt".format(self.seq.description) : self.seq.synced("tcgcgcgtttcggtgatgacggtgaaaacctctg").format("gb").decode("utf-8"),
                           "{}_plan.rst".format(self.seq.description) : plan}
            self.code =  "from {i} import {i} as {c}".format(i=self.seq.description, c="{}")

        else:
            raise Exception("Either one or three arguments!")

    def left_product(self):
        return pydna.pcr(p577, p778, self.seq)

    def assembly_product(self):
        return pydna.pcr(p775, p778, self.seq)

    def right_product(self):
        return pydna.pcr(p775, p578, self.seq)

    def __repr__(self):
        return  "pYPKa_tp_gene_tp({})".format(len(self.seq))





####################################################################################################################


class PathWay(object):

    pYPKpw_lin = pYPKpw.linearize(EcoRV)

    with codecs.open(resource_filename('ypkpathway',os.path.join('data','template_pYPK0_pw.txt')), "r", "utf-8" ) as f:
        codetempl = f.read()

    with codecs.open(resource_filename('ypkpathway',os.path.join('data','template_pYPK0_pw_plan.txt')), "r", "utf-8" ) as f:
        plantempl = f.read()

    def __init__(self, pth, email =  ""):

        self.pth = pth

        self.files = None
        self.tp_gene_tp = []

        pth=self.pth

        while pth:
            last = pth.pop()
            if cloned(pYPK0, (ZraI, EcoRV),  last): # sequence is a tp-gene-tp
                self.tp_gene_tp.append(pYPK0_tp_gene_tp(last))
                continue
            try:
                middle = pth.pop()
                first  = pth.pop()

            except IndexError:
                raise Exception("not enough sequences")
            self.tp_gene_tp.append(pYPK0_tp_gene_tp(first, middle, last))

        self.tp_gene_tp.reverse()

        #self.pYPK0_E_Z = self.tp_gene_tp[0].pYPK0_E_Z
        #pw = [self.pYPK0_E_Z, self.tp_gene_tp[0].left_product()]

        pw = [PathWay.pYPKpw_lin, self.tp_gene_tp[0].left_product()]

        pw.extend([c.assembly_product() for c in self.tp_gene_tp[1:-1]])
        pw.append(self.tp_gene_tp[-1].right_product())

        self.assembly = pydna.Assembly(pw, limit=167-47-10)
        self.seq = self.assembly.circular_products[0]
        self.seq.id = "pYPK0_pathway"
        self.seq.description = self.seq.id

    def theoretical_pathway(self):
        a=pydna.Dseqrecord("GTCgaggaacgccaggttgcccactttctcactagtgacctgcagccGAC")
        b=pydna.Dseqrecord("GTGccatctgtgcagacaaacgcatcagGAT")

        pw=(
        ScTEF1tp,
        a,
        XYL1,
        b,
        ScTDH3tp,
        a,
        XYL2,
        b,
        ScPGI1tp,
        a,
        ScXK,
        b,
        ScFBA1tp,
        a,
        ScTAL1,
        b,
        ScPDC1tp,)

        x=pydna.Dseqrecord("")
        for p in pw:
            x+=p

        return (pCAPs_pSU0_E_Z + x).looped().synced("tcgcgcgtttcggtgatgacggtgaaaacc")

    def generate_files(self):

        self.files = {}

        #self.files["pYPK0_E_Z.txt"] = self.pYPK0_E_Z.format("gb")

        for file_ in (  "voidspace.css",
                        "pYPK0.py",         "pYPK0.txt",
                        "pYPKa.py",         "pYPKa.txt",      "pCAPs.py",
                        "pCAPs.txt",        "pSU0.py",        "pSU0.txt",
                        "pSU0_EcoRV.txt",   "pSU0_EcoRV.py",  "pYPKpw.py",
                        "pYPKpw.txt",       "pYPKpw_lin.txt"                  ):
            with codecs.open( resource_filename('ypkpathway', os.path.join('data', file_)), "r", "utf-8") as f:
                self.files[file_] = f.read()


        date_ = strftime("%Y-%m-%d", gmtime())
        time_ = strftime("%H:%M:%S", gmtime())
        now =  "{} {}".format(date_, time_)

        pathway_name =  "pYPK0_"

        self.files["report.rst"]= ""

        seguids=[]
        pl=  "Specific Primers:\n\n"

        for c in self.tp_gene_tp:
            self.files.update(c.files)
            pathway_name+="{}tp_{}_".format(c.tp1_description.split("tp")[0], c.gene_description)

            self.files["report.rst"]+=  "`{0} <./{0}.txt>`_ (`plan <./{0}_plan.html>`__)\n\n".format(c.seq.description)
            for p in c.pYPKa_clones:
                self.files.update(p.files)
                self.files["report.rst"]+=  "\t * `{0} <./{0}.txt>`_ (`plan <./{0}_plan.html>`__)\n".format(p.name)
                if not p.flag:
                    if not p.insert.seguid() in seguids:
                        pl+="{}\t{}\n".format(p.insert.forward_primer.description, p.insert.forward_primer.seq)
                        pl+="{}\t{}\n".format(p.insert.reverse_primer.description, p.insert.reverse_primer.seq)
                        seguids.append(p.insert.seguid())
            self.files["report.rst"]+= "\n"

        pathway_name+="{}tp_pw".format(self.tp_gene_tp[-1].tp2_description.split("tp")[0])

        self.files["{}.txt".format(pathway_name)] = self.seq.synced("tcgcgcgtttcggtgatgacggtgaaaacctctg").format("gb").decode("utf-8")

        code = []
        cassettes = []

        for i, tgt in enumerate(self.tp_gene_tp):
            cas =  "cas{}".format(i+1)
            code.append(tgt.code.format(cas))
            cassettes.append(cas)

        code.append("cas1  = pcr( p577, p778, cas1)")

        for i in range(2, len(self.tp_gene_tp)):
            code.append("cas{i}  = pcr( p775, p778, cas{i})".format(i=i))
        code.append("cas{i} = pcr( p775, p578, cas{i})".format(i=len(self.tp_gene_tp)))
        code = "\n".join(code)
        self.files["{}.py".format(pathway_name)] = PathWay.codetempl.format(code      = code,
                                                                            cassettes =  ",".join(cassettes),
                                                                            vn = pathway_name )

        pcrs=""

        for frag in self.seq.source_fragments[1:]:
            pcrs +=  "primers {f}, {r} and `{t} <./{t}.txt>`__ => `{n} <./{n}.txt>`__ |br| |br|\n".format(n=frag.name,
                                                                                            f = frag.forward_primer.name,
                                                                                            r = frag.reverse_primer.name,
                                                                                            t = frag.template.description)

            self.files["{n}.txt".format(n=frag.name)] = frag.format("gb").decode("utf-8")

        self.files["{}_plan.rst".format(pathway_name)] = PathWay.plantempl.format( line =  "="*len(pathway_name),
                                                                                   name= pathway_name,
                                                                                   number = len(self.tp_gene_tp),
                                                                                   figure  = add_space(self.seq.small_fig()),
                                                                                   pcrs = pcrs)

        with codecs.open(resource_filename('ypkpathway',os.path.join('data','template_primer_list.txt')), "r", "utf-8" ) as f:
            self.files["primer_list.txt"] = f.read().format(primer_list = pl)

        with codecs.open(resource_filename('ypkpathway',os.path.join('data','template_report.txt')), "r", "utf-8" ) as f:
            self.files["report.rst"] = f.read().format(**locals())+self.files["report.rst"]

        #print type(self.files["report.rst"]);import sys;sys.exit()

        args = {'output_encoding' : 'unicode',
                'input_encoding'  : 'unicode',
                'stylesheet_path' : resource_filename('ypkpathway', os.path.join('data','voidspace.css'))}

        for f, c in self.files.items():
            if f.endswith(".rst"):
                html = publish_string(c,
                                      writer=HisWriter(),
                                      settings=None,
                                      settings_overrides=args)

                self.files[f.split(".")[0]+".html"] = html


    def zs(self):
        if not self.files:
            self.generate_files()
        zipstream = cStringIO.StringIO()
        with zipfile.ZipFile(zipstream, mode='a') as myzip:
            for name, content in self.files.items():
                myzip.writestr(name, content)
        zipstream.seek(0)
        return zipstream.getvalue()

def main():
    import docopt
    try:
        arguments = docopt.docopt(__doc__)
    except docopt.DocoptExit as e:
        print e.message

    if arguments['--version']:
        from ._version import get_versions
        __version__      = get_versions()['version'][:5]
        del get_versions
        print "ypkpathway version:",__version__

    if arguments['<path>']:
        file_ = arguments['<path>']
        try:
            with open(file_, "rU") as f: text=f.read()
        except IOError:
            print arguments['<path>'], 'could not be opened!'
            sys.exit(1)

        print "Assembly started! (This might take a while...)"

        pw = PathWay( pydna.parse( text ) )
        pw.generate_files()

        try:
            os.makedirs('ypk_assembly')
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise

        for name, content in pw.files.items():
            with codecs.open(os.path.join('ypk_assembly',name),'w', "utf8") as f:
                f.write(content)
        print "Assembly finished! files written to folder ypkpathway"

if __name__ == "__main__":

    cwd = os.getcwd()
    abspath = os.path.abspath(__file__)
    dname = os.path.dirname(abspath)
    os.chdir(dname)
    main()
    os.chdir(cwd)

