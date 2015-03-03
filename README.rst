==========
ypkpathway
==========

.. image:: https://travis-ci.org/BjornFJohansson/ypkpathway.svg 
    :target: https://travis-ci.org/BjornFJohansson/ypkpathway
    
.. image:: https://coveralls.io/repos/BjornFJohansson/ypkpathway/badge.svg?branch=master 
    :target: https://coveralls.io/r/BjornFJohansson/ypkpathway?branch=master
  
.. image:: https://readthedocs.org/projects/ypkpathway/badge/?version=latest
    :target: https://readthedocs.org/projects/ypkpathway/?badge=latest
    :alt: Documentation Status

.. image:: https://pypip.in/download/ypkpathway/badge.svg
    :target: https://pypi.python.org/pypi/ypkpathway/
    :alt: Downloads
    
.. image:: https://pypip.in/version/ypkpathway/badge.svg
    :target: https://pypi.python.org/ypkpathway/
    :alt: Latest Version

.. image:: https://pypip.in/wheel/ypkpathway/badge.svg
    :target: https://pypi.python.org/pypi/ypkpathway/
    :alt: Wheel Status
    
Ypkpatwhay provides a command line application for planning DNA assembly projects 
using the Yeast Pathway Kit protocol. 

Typical usage at the command line could look like this::

    bjorn@bjorn-UL30A:/$ ypkpathway four_gene_xylose_pathway.txt
    
Where four_gene_xylose_pathway1.txt is a text file containing DNA sequences to be assembled
in FASTA or Genbank format. See documentation for file format.

The ypkpathway creates a subdirectory called generates a sub directory called "ypk_assembly".
This subdirectory will contain a file summary of the DNA assembly project called report.html 
which can be opened in a web browser. 

The report will include:

* The sequence of the final pathway
* Sequences of all generated intermediate vectors
* All PCR primers needed for the amplification of pathway components
* Diagnostic PCR product fragment lengths indicating correct and incorrect clonings


ypkpathway depends on the `pydna <https://pypi.python.org/pypi/pydna/>`_ package.



NEWS
====

=======   ========== =============================================================
version   date       comment
=======   ========== =============================================================
0.5.0	  2014-12-09 Changed to work with new caching version of pydna. Fixed 
				     unicode errors.

0.0.9     2014-06-14 bugfix
                     changed errors in automatically generated pydna source files
                     
0.0.8     2014-05-14 first release
=======   ========== =============================================================

System Requirements
===================

- `Python 2.7 <http://www.python.org>`_.

- `pydna>=0.7.2 <https://pypi.python.org/pypi/pydna/>`_.

- `docutils>=0.11 <https://pypi.python.org/pypi/docutils/>`_.

- `docopt>=0.6.1 <https://pypi.python.org/pypi/docopt/>`_.


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

Make sure all dependencies are intalled. Open the pydna source code 
directory (containing the setup.py file) in terminal and type::

    sudo python setup.py install <enter>

If you need to do additional configuration, e.g. changing the base
directory, please type `python setup.py`, or see the documentation for
Setuptools.

Distribution Structure
======================

README.txt          -- This file.

LICENSE.txt         -- What you can do with the code.

MANIFEST.in         -- Tells distutils what files to distribute

setup.py            -- Installation file.

run_tests.py        -- run tests by "python run_tests.py"<enter> (Currently empty)

ypkpathway/         -- The actual code.

docs/               -- Documentation.

tests/              -- Testing code and data.
