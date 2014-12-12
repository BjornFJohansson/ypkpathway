#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pydna import Genbank
gb = Genbank("bjornjobbb@gmail.com")
pCAPs = gb.nucleotide("AJ001614")
pCAPs.features = pCAPs.features [1:]

