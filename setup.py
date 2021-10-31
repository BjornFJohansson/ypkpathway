#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""docstring."""

from setuptools import setup
from setuptools import Command
from setuptools import find_packages
from os import path
import re

# Read __author__, __email__. from __init__.py
__author__ = "__author__"
__email__ = "__email__"
for line in open('src/ypkpathway/__init__.py'):
    if line.startswith('__') and not line.startswith('__version') and not line.startswith('__long'):
        exec(line.strip())

this_directory = path.abspath(path.dirname(__file__))
with open(path.join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


install_requires = []
with open("requirements.txt", encoding="utf-8") as f:
    for line in f.readlines():
        install_requires.append(re.split(r"(<|>|=)=", line)[0])


class PyTest(Command):
    """doctsring."""

    user_options = []

    def initialize_options(self):
        """doctsring."""
        pass

    def finalize_options(self):
        """doctsring."""
        pass

    def run(self):
        """doctsring."""
        import subprocess
        import sys

        errno = subprocess.call([sys.executable, "run_test.py"])
        raise SystemExit(errno)


setup(  name='ypkpathway',
        author          = __author__,
        author_email    = __email__,
        zip_safe = False,
        cmdclass={"test": PyTest},
        packages=find_packages("src"),
        package_dir={"": "src"},
        url='https://github.com/BjornFJohansson/ypkpathway',
        license='LICENSE.txt',
        description='Simulation and documentation of metabolic pathway assemblies using the Yeast Pathway Kit.',
        long_description=long_description,
        long_description_content_type='text/markdown',
        setup_requires=["pytest-runner", "setuptools_scm"],
		tests_require=["pytest"],
		use_scm_version={"write_to": "src/ypkpathway/_version.py"},
		install_requires=install_requires,
		keywords="bioinformatics",
        package_data={'ypkpathway': [path.join('data','*'), path.join('icons','*')]},
        entry_points = { 'gui_scripts'     : ['ypkpathwaygui = ypkpathway.gui:main',],
                         'console_scripts' : ['ypkpathway    = ypkpathway.ypkpathway:main']},
        classifiers = ['Development Status :: 4 - Beta',
                       'Environment :: Console',
                       'Intended Audience :: Education',
                       'Intended Audience :: Science/Research',
                       'License :: OSI Approved :: BSD License',
                       'Programming Language :: Python :: 3.7',
                       'Programming Language :: Python :: 3.8',
                       'Programming Language :: Python :: 3.9',
                       'Topic :: Education',
                       'Topic :: Scientific/Engineering :: Bio-Informatics',])
