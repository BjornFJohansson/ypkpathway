==========
ypkpathway
==========

.. image:: https://travis-ci.org/BjornFJohansson/ypkpathway.svg
    :target: https://travis-ci.org/BjornFJohansson/ypkpathway

.. image:: https://ci.appveyor.com/api/projects/status/ol5ru8po7fx5cimj?svg=true
    :target: https://ci.appveyor.com/project/BjornFJohansson/ypkpathway

.. image:: https://img.shields.io/pypi/v/ypkpathway.png
    :target: https://pypi.python.org/pypi/ypkpathway/
    :alt: Downloads

.. image:: https://img.shields.io/pypi/dm/ypkpathway.png
    :target: https://pypi.python.org/pypi/ypkpathway/
    :alt: Latest Version

.. image:: https://www.versioneye.com/user/projects/55645b646361300021ae0200/badge.svg?style=flat(Dependency Status)!
    :target: https://www.versioneye.com/user/projects/55645b646361300021ae0200
    :alt: versioneye


Ypkpatwhay provides a command line application for planning DNA assembly projects
using the Yeast Pathway Kit protocol. 

Typical usage at the command line could look like this::

    bjorn@bjorn-UL30A:/$ ypkpathway four_gene_xylose_pathway.txt

Where four_gene_xylose_pathway1.txt is a text file containing DNA sequences in FASTA or Genbank format 
to be assembled. The ypkpathway command above creates a folder with a series of
IPython notebooks describing the assembly process simulated with pydna. Help is available by the -h option:

    bjorn@bjorn-UL30A:/$ ypkpathway -h
    Usage: ypkpathway <path> [<dir>]
           ypkpathway -h|--help
           ypkpathway -v|--version
           ypkpathway -t|--test

    Arguments:
        <path>  path to data file containing sequences to be assembled

        <dir>   Directory to put generated sequence files,defaults to
                <ypk_assembly> in the current working directory.

    Options:
        -h, --help      Show this screen.
        -v, --version   Show version.
        -t, --test      Run tests.


The results folder will include:

* The sequence of the final pathway and all intermediate vectors in Genbank format
* IPython notebooks file describing the final assembly and intermediate assemblies.
* All PCR primers needed for the amplification of pathway components.
* Expected diagnostic PCR product fragment lengths indicating correct and incorrect clonings.


NEWS
====

=======   ========== =============================================================
version   date       comment
=======   ========== =============================================================
0.8.4     2015-05-30 An IPython.display.FileLink object is now returned from the 
                     ypkpathway function.

0.8.3     2015-05-29 windows fix for opening final pathway notebook.

0.8.2     2015-05-28 pydna dependency 0.9.2, fix a unicode env variable.

0.8.1     2015-05-26 pydna dependency 0.9.1

0.8.0     2015-05-26 Complete change in functionality. IPython notebook files are
                     generated and executed instead of rst and html files.

0.6.0     2015-04-17 Added nosetests

0.5.0	  2014-12-09 Changed to work with new caching version of pydna. Fixed
				     unicode errors.

0.0.9     2014-06-14 bugfix
                     changed errors in automatically generated pydna source files

0.0.8     2014-05-14 first release
=======   ========== =============================================================

System Requirements
===================

- `Python 2.7 <http://www.python.org>`_.

- `pydna>=0.9.1 <https://pypi.python.org/pypi/pydna/>`_.

- `ipython>=3.1.0 <https://pypi.python.org/pypi/ipython/>`_.

- `docopt>=0.6.2 <https://pypi.python.org/pypi/docopt/>`_.

- `notedown>=1.4.4 <https://pypi.python.org/pypi/notedown/>`_.


Python 2.x
----------

Versions other than 2.7 has not been tried with this software.
Version 2.7.3 was used to build the distribution.

Python 3.x
----------

This code has not been tested with python 3.

Installation
============

The best way of installing ypkpathway is with `pip <https://pypi.python.org/pypi/pip/>`_::

    sudo pip install ypkpathway <enter>

Pip will take care of the installation of any missing dependencies.

Source
------

Make sure all dependencies are installed. Open the pydna source code
directory (containing the setup.py file) in terminal and type::

    sudo python setup.py install <enter>

If you need to do additional configuration, e.g. changing the base
directory, please type `python setup.py`, or see the documentation for
Setuptools.

Testing
-------

Nose is needed to run the test suit.

    sudo pip install nose <enter>

Then from the source directory do:

    python run_tests.py

Distribution Structure
======================

README.rst          -- This file.

LICENSE.txt         -- What you can do with the code.

MANIFEST.in         -- Tells distutils what files to distribute

setup.py            -- Installation file.

run_tests.py        -- run tests by "python run_tests.py"<enter>

ypkpathway/         -- The actual code.

docs/               -- Documentation.

tests/              -- Testing code and data.
