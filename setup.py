#!/usr/bin/env python
# -*- coding: utf-8 -*-

import versioneer
versioneer.VCS = 'git'
versioneer.versionfile_source = 'ypkpathway/_version.py'
versioneer.versionfile_build = 'ypkpathway/_version.py'
versioneer.tag_prefix = '' # tags are like 1.2.0
versioneer.parentdir_prefix = '' # dirname like 'myproject-1.2.0'

# Read author etc..
for line in open('ypkpathway/__init__.py'):
    if line.startswith('__') and not line.startswith('__version') and not line.startswith('__long'):
        exec(line.strip())

from setuptools import setup
import os

setup(  name='ypkpathway',
        version=versioneer.get_version()[:5],
        cmdclass=versioneer.get_cmdclass(),
        author          =__author__,
        author_email    =__email__,
        packages=['ypkpathway'],
        package_data={'ypkpathway': [os.path.join('data','*')]},
        entry_points = { 'console_scripts': [ 'ypkpathway = ypkpathway.ypkpathway:main' ]},
        url='http://pypi.python.org/pypi/ypkpathway/',
        license='LICENSE.txt',
        description='''Assemble metabolic pathways from the command line using the ypkpathway protocol''',
        long_description=open('README.rst').read(),
        install_requires =[ "pydna", "ipython", "docopt", "notedown"],
        test_suite="run_tests.load_my_tests",
        zip_safe = False,
        keywords = "bioinformatics",
        classifiers = ['Development Status :: 4 - Beta',
                       'Environment :: Console',
                       'Intended Audience :: Education',
                       'Intended Audience :: Science/Research',
                       'License :: OSI Approved :: BSD License',
                       'Programming Language :: Python :: 2.7',
                       'Topic :: Education',
                       'Topic :: Scientific/Engineering :: Bio-Informatics',])
