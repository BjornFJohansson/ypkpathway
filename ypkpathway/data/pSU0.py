#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydna import Genbank
gb = Genbank("bjornjobbb@gmail.com")
pSU0 = gb.nucleotide("AB215109")
pSU0.features = pSU0.features[1:]
