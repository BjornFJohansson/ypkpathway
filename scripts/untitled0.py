'''
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
