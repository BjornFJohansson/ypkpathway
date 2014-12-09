try:
    seq = read(fname)
except ValueError:
    print "Could not read {}".format(fname)
    print " Executing     {}".format(fname)
    globals().update({name:design()})
    globals()[name].write(fname)
    print "Written        {}".format(fname)
    seq = None
else:
    print "Reading        {}".format(fname)
    globals().update({name:seq})

if __name__=="__main__" or "redo" in sys.argv:
    from pydna import eq
    import shutil
    #print " Executing     {}".format(fname) 
    if seq:
        new = design()
        if str(new.seq).lower()==str(seq.seq).lower():
            print "new and old sequences are identical."
        else:
            print "new ({}bp) and old ({}bp) sequences are NOT identical!".format(len(new),len(seq))
            if eq(new, seq):
                print "new and old sequences are equivalent."
            else:
                print "new ({}bp) and old ({}bp) sequences are NOT equivalent!".format(len(new),len(seq))
            oldname = "{}_old.gb".format(name)
            shutil.move(fname, oldname)
            new.write(fname)
            print "old = {}".format(oldname)
            print "new = {}".format(fname)
