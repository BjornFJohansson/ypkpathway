#!/usr/bin/env python
# -*- coding: utf-8 -*-

import versioneer
#versioneer.VCS = 'git'
#versioneer.versionfile_source = 'ypkpathway/_version.py'
#versioneer.versionfile_build = 'ypkpathway/_version.py'
#versioneer.tag_prefix = '' # tags are like 1.2.0
#versioneer.parentdir_prefix = '' # dirname like 'myproject-1.2.0'

# Read author etc..
for line in open('ypkpathway/__init__.py'):
    if line.startswith('__') and not line.startswith('__version') and not line.startswith('__long'):
        exec(line.strip())

from setuptools import setup
import os#, codecs

setup(  name='ypkpathway',
        version=versioneer.get_version()[:5],
        cmdclass=versioneer.get_cmdclass(),
        author          =__author__,
        author_email    =__email__,
        packages=['ypkpathway'],
        package_data={'ypkpathway': [os.path.join('data','*'), os.path.join('icons','*')]},
        entry_points = { 'gui_scripts'    : [ 'ypkpathway = ypkpathway.gui:main',   ], 'console_scripts' : ['ypkpathway_cli = ypkpathway.ypkpathway:main']    },
        url='http://pypi.python.org/pypi/ypkpathway/',
        license='LICENSE.txt',
        description='''Simulation and documentation of metabolic pathway assembly using the Yeast Pathway Kit.''',
        #long_description=codecs.open('README.rst', "r", "utf8").read(),
        long_description=open('README.rst').read(),
        #setup_requires=['setuptools-markdown'],
        #long_description_markdown_filename='README.md',
        install_requires =[ "pydna", "ipython", "docopt", "notedown"],
        test_suite="run_tests.load_my_tests",
        zip_safe = False,
        keywords = u"bioinformatics",
        classifiers = [u'Development Status :: 4 - Beta',
                       u'Environment :: Console',
                       u'Intended Audience :: Education',
                       u'Intended Audience :: Science/Research',
                       u'License :: OSI Approved :: BSD License',
                       u'Programming Language :: Python :: 2.7',
                       u'Topic :: Education',
                       u'Topic :: Scientific/Engineering :: Bio-Informatics',])







