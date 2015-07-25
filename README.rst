|icon| ypkpathway
=================

Ypkpathway is a simulator and documentation generator for *in-vivo*
pathway assembly using the Yeast Pathway Kit protocol. It takes as
argument a series of sequences in a text file and produces a self
contained folder containing assembled sequences of intermediary vectors,
final assembly and PCR primers. Other useful information such as PCR
conditions and PCR product sizes are also included. The documantation of
the assembly is given in the narrative IPython notebook format which can
be executed independently of ypkpathway.

The assembly process is simulated using
`pydna <https://github.com/BjornFJohansson/pydna>`__ in `IPython
notebooks <http://ipython.org/notebook.html>`__ which are automatically
generated and executed.

See an example of an four gene assembly
`here <http://nbviewer.ipython.org/github/BjornFJohansson/ypkpathway/blob/master/docs/pth6/pw.ipynb>`__.
The notebooks in the example are located in the docs folder in this
repository and visualized through
`nbviewer <http://nbviewer.ipython.org/>`__. The example above was made
with the
`pth6 <http://nbviewer.ipython.org/github/BjornFJohansson/ypkpathway/blob/master/tests/pth6.txt>`__
indata.

There are five more example indata example files that are a part of the
automatic test suit:
`pth1 <http://nbviewer.ipython.org/github/BjornFJohansson/ypkpathway/blob/master/tests/pth1.txt>`__
`pth2 <http://nbviewer.ipython.org/github/BjornFJohansson/ypkpathway/blob/master/tests/pth2.txt>`__
`pth3 <http://nbviewer.ipython.org/github/BjornFJohansson/ypkpathway/blob/master/tests/pth3.txt>`__
`pth4 <http://nbviewer.ipython.org/github/BjornFJohansson/ypkpathway/blob/master/tests/pth4.txt>`__
`pth5 <http://nbviewer.ipython.org/github/BjornFJohansson/ypkpathway/blob/master/tests/pth5.txt>`__
`pth7 <http://nbviewer.ipython.org/github/BjornFJohansson/ypkpathway/blob/master/tests/pth7.txt>`__

|screenshot|

|pypi| Python wheel and source distributions on Pypi.

|travis| Python wheel and source distributions on Pypi are built and
tested on on travis-ci using a linux back end.

|appveyor| 64 bit Windows binary executable and a Windows binary wheel
are built on Appveyor-CI

|anaconda build| Anaconda packages for Windows, Linux and MacOSX are
built at anaconda.org on a Linux back end.

|dependencies| Dependencies are tracked at versioneye.

|pypi dl| Pypi download count.

|anaconda| Anaconda download count.

|GitHub tag| Github repository.

The ypkpatwhay package provides a graphical point and click interface
and a command line application for planning DNA assembly projects using
the Yeast Pathway Kit protocol.

Please refer to the
`manual <https://github.com/BjornFJohansson/ypkpathway/blob/master/docs/manual.pdf>`__
for details on how to use the software.

Installation
------------

The best way of installing ypkpathway is by first installing the free
`Anaconda Python
distribution <https://store.continuum.io/cshop/anaconda/>`__ which comes
with many packages and dependencies out of the box. Using the conda
package manages simply type:

::

    bjorn@bjorn-UL30A:/$ conda install ypkpathway

and the app and all dependencies will be installed. The
`manual <https://github.com/BjornFJohansson/ypkpathway/blob/master/docs/manual.pdf>`__
contains a detailed walk through of this installation option.

Alternatively, ypkpathway can be installed using
`pip <https://pypi.python.org/pypi/pip>`__ which is the `PyPA
recommended <https://python-packaging-user-guide.readthedocs.org/en/latest/current.html>`__
tool for installing Python packages.

::

    bjorn@bjorn-UL30A:/$ pip install ypkpathway

Pip may have trouble to install two dependecies
`biopython <https://pypi.python.org/pypi/biopython>`__ which is a
dependency of pydna and
`PyQt4 <https://pypi.python.org/pypi/PyQt4/4.11.4>`__ which have binary
extensions. These can be separately installed. Binary installers of
PyQt4 can be found
`here <http://www.riverbankcomputing.com/software/pyqt/download>`__.
Instructions for how to install Biopython can be found
`here <http://biopython.org/wiki/Download>`__.

Ypkpathway can also be installed from source by downloading one of the
source distributions. Unpack the zip or .tar.gz archive and type:

::

    bjorn@bjorn-UL30A:/$ python setup.py install

Dependencies has to be manually installed in this case. There are also
.exe installers for ypkpathway that can be installed by double clicking.
These do not install the dependencies either.

Dependencies
------------

The ypkpathway dependencies are pure Python modules except for PyQt.
Pydna depends on `biopython <https://pypi.python.org/pypi/biopython>`__
which has to be installed using a binary installer or a c-compiler has
to be present.

`pydna <https://pypi.python.org/pypi/pydna>`__

`ipython <https://pypi.python.org/pypi/ipython>`__

`docopt <https://pypi.python.org/pypi/docopt>`__

`notedown <https://pypi.python.org/pypi/notedown>`__

`PyQt4 <https://pypi.python.org/pypi/PyQt4>`__

`tornado <https://pypi.python.org/pypi/tornado>`__

`terminado <https://pypi.python.org/pypi/terminado>`__

`pyzmq <https://pypi.python.org/pypi/pyzmq>`__

There are binary Windows installers for IPython, Bioppython, docopt
`here <http://www.lfd.uci.edu/~gohlke/pythonlibs/>`__

Graphical user interface
------------------------

The application is written in PyQt4 and can be started from the command
line by typing ypkpathway and pressing :

::

    bjorn@bjorn-UL30A:/$ ypkpathway

It can also be started from the Anaconda Launcher if installed using
conda on the Anaconda Python distribution.

Command line interface
----------------------

Typical usage at the command line could look like this:

::

    bjorn@bjorn-UL30A:/$ ypkpathway_cli pth6.txt

Where pth6.txt is a text file containing DNA sequences in FASTA or
Genbank format to be assembled as described in the
`manual <https://github.com/BjornFJohansson/ypkpathway/blob/master/docs/manual.pdf>`__.

The ypkpathway\_cli command above creates a folder with a series of
IPython notebooks describing the assembly process simulated with pydna.
Help is available by the -h option:

::

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

Results
-------

The ypkpathway and ypkpathway\_cli both produce the same result, which
is a results folder containing a selection of files. The folder will
contain:

-  The sequence of the final pathway and all intermediate vectors in
   `Genbank <http://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html>`__
   format
-  IPython notebooks files describing the final assembly and
   intermediate assemblies.
-  All PCR primers needed for the amplification of pathway components.
-  Expected diagnostic PCR product fragment lengths indicating correct
   and incorrect clonings.

The IPython notebook files in the results folder can be viewed with a
web browser with oly IPython is installed on the computer. There are
static versions of the notebook files that can be viewed with only a web
browser (not eve Python is required).

Development
-----------

Ypkpathway is open source software and developen on Github |GitHub tag|.

.. |icon| image:: https://raw.githubusercontent.com/BjornFJohansson/ypkpathway/master/icon.resized.png
   :target: https://pypi.python.org/pypi/ypkpathway/
.. |screenshot| image:: https://raw.githubusercontent.com/BjornFJohansson/ypkpathway/master/screenshot.resized.png
   :target: https://github.com/BjornFJohansson/ypkpathway
.. |pypi| image:: https://img.shields.io/pypi/v/ypkpathway.png
   :target: https://pypi.python.org/pypi/ypkpathway/
.. |travis| image:: https://travis-ci.org/BjornFJohansson/ypkpathway.svg
   :target: https://travis-ci.org/BjornFJohansson/ypkpathway
.. |appveyor| image:: https://ci.appveyor.com/api/projects/status/ol5ru8po7fx5cimj?svg=true
   :target: https://ci.appveyor.com/project/BjornFJohansson/ypkpathway
.. |anaconda build| image:: https://anaconda.org/bjornfjohansson/ypkpathway/badges/build.svg
   :target: https://anaconda.org/bjornfjohansson/ypkpathway/builds
.. |dependencies| image:: https://www.versioneye.com/user/projects/55645b646361300021ae0200/badge.svg?style=flat(Dependency%20Status)!
   :target: https://www.versioneye.com/user/projects/55645b646361300021ae0200
.. |pypi dl| image:: https://img.shields.io/pypi/dm/ypkpathway.png
   :target: https://pypi.python.org/pypi/ypkpathway/
.. |anaconda| image:: https://anaconda.org/bjornfjohansson/ypkpathway/badges/downloads.svg
   :target: https://anaconda.org/bjornfjohansson/ypkpathway
.. |GitHub tag| image:: https://img.shields.io/github/tag/BjornFJohansson/ypkpathway.svg
   :target: https://github.com/BjornFJohansson/ypkpathway
