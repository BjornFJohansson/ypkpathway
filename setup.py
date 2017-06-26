#!/usr/bin/env python
# -*- coding: utf-8 -*-

import versioneer

# Read author etc. from __init__.py
for line in open('ypkpathway/__init__.py'):
    if line.startswith('__') and not line.startswith('__version') and not line.startswith('__long'):
        exec(line.strip())

from setuptools import setup

import os

with open("README.md") as f: long_description = f.read()
    
try:
    from pypandoc import convert_file
except ImportError:
    print("warning: pypandoc module not found, could not convert Markdown to RST")
    with open("README.md", encoding="utf-8") as f:
        long_description = f.read()
else:
    long_description = "\n"+convert_file("README.md", 'rst')   

setup(  name='ypkpathway',
        version=versioneer.get_version()[:5],
        cmdclass=versioneer.get_cmdclass(),
        author          =__author__,
        author_email    =__email__,
        packages=['ypkpathway'],
        package_data={'ypkpathway': [os.path.join('data','*'), os.path.join('icons','*')]},
        entry_points = { 'gui_scripts'     : [ 'ypkpathway = ypkpathway.gui:main',   ], 
                         'console_scripts' : [ 'ypkpathway_cli = ypkpathway.ypkpathway:main']    },
        url='http://pypi.python.org/pypi/ypkpathway/',
        license='LICENSE.txt',
        description='''Simulation and documentation of metabolic pathway assemblies using the Yeast Pathway Kit.''',
        long_description=long_description,
        
        setup_requires=['pytest-runner', 'versioneer', 'pypandoc'],
        tests_require =['pytest', "pydna", "ipython", "docopt", "notedown", "nbformat", "nbconvert", "pydna", "notedown"],
        install_requires = [ "pydna", "ipython", "docopt", "notedown", "nbformat", "nbconvert", "pydna", "notedown"],

        zip_safe = False,
        keywords = u"bioinformatics",
        classifiers = ['Development Status :: 4 - Beta',
                       'Environment :: Console',
                       'Intended Audience :: Education',
                       'Intended Audience :: Science/Research',
                       'License :: OSI Approved :: BSD License',
                       'Programming Language :: Python :: 3.5',
                       'Programming Language :: Python :: 3.6',
                       'Topic :: Education',
                       'Topic :: Scientific/Engineering :: Bio-Informatics',])
