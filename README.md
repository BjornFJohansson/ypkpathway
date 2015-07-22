#[![icon](https://raw.githubusercontent.com/BjornFJohansson/ypkpathway/master/icon.resized.png)](https://pypi.python.org/pypi/ypkpathway/) ypkpathway 


Ypkpathway is a simulator and documentation generator for _in-vivo_ pathway assembly using the Yeast Pathway Kit protocol. It takes as argument a series of
sequences in a text file and produces a self contained folder containing assembled sequences of intermediary vectors, final assembly and PCR primers. 
Other useful information such as PCR conditions and PCR product sizes are also included. The documantation of the assembly is given in the narrative 
IPython notebook format which can be executed independently of ypkpathway.

The assembly process is simulated using [pydna](https://github.com/BjornFJohansson/pydna) in [IPython notebooks](http://ipython.org/notebook.html) which are automatically generated and executed. 

[![screenshot](https://raw.githubusercontent.com/BjornFJohansson/ypkpathway/master/screenshot.resized.png)](https://github.com/BjornFJohansson/ypkpathway)

A wheel and setuptools pypi [![pypi](https://img.shields.io/pypi/v/ypkpathway.png)](https://pypi.python.org/pypi/ypkpathway/)
 source distributions are built and tested on on travis-ci using a linux back end [![travis](https://travis-ci.org/BjornFJohansson/ypkpathway.svg)](https://travis-ci.org/BjornFJohansson/ypkpathway). A 64 bit Windows binary executable and a Windows binary wheel are built on Appveyor CI 
[![appveyor](https://ci.appveyor.com/api/projects/status/ol5ru8po7fx5cimj?svg=true)](https://ci.appveyor.com/project/BjornFJohansson/ypkpathway) on a Windows back end. Binstar builds [![anaconda build](https://anaconda.org/bjornfjohansson/ypkpathway/badges/build.svg)](https://anaconda.org/bjornfjohansson/ypkpathway/builds) Anaconda packages for Windows, Linux and MacOSX on a Linux back end [![anaconda](https://anaconda.org/bjornfjohansson/ypkpathway/badges/downloads.svg)](https://anaconda.org/bjornfjohansson/ypkpathway).

Pypi downloads [![pypi dl](https://img.shields.io/pypi/dm/ypkpathway.png)](https://pypi.python.org/pypi/ypkpathway/). Dependencies are tracked at versioneye 
[![dependencies](https://www.versioneye.com/user/projects/55645b646361300021ae0200/badge.svg?style=flat(Dependency%20Status)!)](https://www.versioneye.com/user/projects/55645b646361300021ae0200).


The ypkpatwhay package provides a graphical point and click interface and a command line application for planning DNA assembly projects 
using the Yeast Pathway Kit protocol.


Please refer to the [manual](https://github.com/BjornFJohansson/ypkpathway/blob/master/docs/manual.pdf) for details on how to use the software.


##Installation

The best way of installing ypkpathway is by first installing the free [Anaconda distribution](https://store.continuum.io/cshop/anaconda/) which comes with
many packages and dependencies out of the box. Using the conda package manages simply type:

    bjorn@bjorn-UL30A:/$ conda install ypkpathway

and the app and all dependencies will be installed. The [manual](https://github.com/BjornFJohansson/ypkpathway/blob/master/docs/manual.pdf) contains a detailed 
walk through of this installation option.

Alternatively, ypkpathway can be installed using [pip](https://pypi.python.org/pypi/pip) which is the [PyPA recommended](https://python-packaging-user-guide.readthedocs.org/en/latest/current.html) tool for installing Python packages.

    bjorn@bjorn-UL30A:/$ pip install ypkpathway

Pip may have trouble to install two dependecies [biopython](https://pypi.python.org/pypi/biopython/1.65) which is a dependency of pydna and [PyQt4](https://pypi.python.org/pypi/PyQt4/4.11.4) which have binary extensions. 
These can be separately installed. Binary installers of PyQt4 can be found [here](http://www.riverbankcomputing.com/software/pyqt/download). Instructions for how to install Biopython can 
be found [here](http://biopython.org/wiki/Download).

Ypkpathway can also be installed from source by downloading one of the source distributions. Unpack the zip or .tar.gz archive and type:

    bjorn@bjorn-UL30A:/$ python setup.py install

Dependencies has to be manually installed in this case. There are also .exe installers for ypkpathway that can be installed by double clicking.
These do not install the dependencies either.

##Dependencies

[pydna](https://pypi.python.org/pypi/pydna)

[ipython](https://pypi.python.org/pypi/ipython)

[docopt](https://pypi.python.org/pypi/docopt)

[notedown](https://pypi.python.org/pypi/notedown)

[PyQt4](https://pypi.python.org/pypi/PyQt4)


There are binary Windows installers for IPython, Bioppython, docopt [here](http://www.lfd.uci.edu/~gohlke/pythonlibs/) 


tornado>=4.1
terminado>=0.5
pyzmq>=14.6.0


##Graphical user interface

The application is written in PyQt4 and can be started from the command line by typing ypkpathway and pressing <enter>:

    bjorn@bjorn-UL30A:/$ ypkpathway

It can also be started from the Anaconda Launcher if installed using conda on the Anaconda Python distribution.


##Command line interface

Typical usage at the command line could look like this:

    bjorn@bjorn-UL30A:/$ ypkpathway_cli pth6.txt

Where pth6.txt is a text file containing DNA sequences in FASTA or Genbank format to be assembled as described in the [manual](https://github.com/BjornFJohansson/ypkpathway/blob/master/docs/manual.pdf).

The ypkpathway_cli command above creates a folder with a series of IPython notebooks describing 
the assembly process simulated with pydna. Help is available by the -h option:

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


##Results

The ypkpathway and ypkpathway_cli both produce the same result, which is a results folder containing a selection of files.
The folder will contain:


-   The sequence of the final pathway and all intermediate vectors in Genbank format
-   IPython notebooks file describing the final assembly and intermediate assemblies.
-   All PCR primers needed for the amplification of pathway components.
-   Expected diagnostic PCR product fragment lengths indicating correct and incorrect clonings.

The IPython notebook files in the results folder can be viewed with a web browser with oly IPython is installed on the computer.
There are static versions of the notebook files that can be viewed with only a web browser (not eve Python is required).

##Development

Ypkpathway is open source software and developen on 
Github [![GitHub tag](https://img.shields.io/github/tag/BjornFJohansson/ypkpathway.svg)](https://github.com/BjornFJohansson/ypkpathway).



