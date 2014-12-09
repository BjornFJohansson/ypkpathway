#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Read version numbers, author etc..
__version__ = "Undefined"
for line in open('./ypkpathway/__init__.py'):
    if line.startswith('__'):
        exec(line.strip())  

from setuptools import setup
import os

setup(  name='ypkpathway',
        version         =__version__,
        author          =__author__,
        author_email    =__email__,
        packages=['ypkpathway'],
        
        package_data={'ypkpathway': [os.path.join('data','*')]},
        
        #package_data={'ypkpathway': [os.path.join('ypkpathway','data','*')]},
        
        entry_points = { 'console_scripts': [ 'ypkpathway = ypkpathway.ypkpathway:main' ]},
        
        url='http://pypi.python.org/pypi/py-genome/',
        license='LICENSE.txt',
        description='''Assemble metabolic pathways from the command line''',
        long_description=open('README.txt').read(),
        install_requires =[ "pydna>=0.6.1", "docutils>=0.11", "docopt>=0.6.1" ],
        test_suite="run_tests.load_my_tests",
        zip_safe = False,
        keywords = "bioinformatics",
        classifiers = ['Development Status :: 3 - Alpha',
                       'Environment :: Console',
                       'Intended Audience :: Education',
                       'Intended Audience :: Science/Research',
                       'License :: OSI Approved :: BSD License',
                       'Programming Language :: Python :: 2.7',
                       'Topic :: Education',
                       'Topic :: Scientific/Engineering :: Bio-Informatics',])
