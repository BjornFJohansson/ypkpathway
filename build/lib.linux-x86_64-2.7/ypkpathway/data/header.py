#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
sys.dont_write_bytecode = True
import os
from pydna import read
from pydna import Dseqrecord
name = os.path.split(os.path.splitext(__file__)[0])[-1]
fname = name+".txt"
