#!/usr/bin/env python
# -*- coding: utf-8 -*-
__author__      ="Björn Johansson"
__date__        ="April 17, 2015"
__copyright__   ="Copyright 2013, 2014, 2015 Björn Johansson"
__credits__     = ["Björn Johansson"]
__license__     = "BSD"
__maintainer__  = "Björn Johansson"
__email__       = "bjorn_johansson@bio.uminho.pt"

from ._version import get_versions
__version__      = get_versions()['version'][:5]
__long_version__ = get_versions()['version']
del get_versions

from .ypkpathway import pathway

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
