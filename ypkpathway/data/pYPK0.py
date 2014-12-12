#!/usr/bin/env python
# -*- coding: utf-8 -*-

from pCAPs import pCAPs
from pSU0_EcoRV import pSU0_EcoRV
from pydna import Assembly2
from Bio.Restriction import BamHI, EcoRI, PvuI

psu = pSU0_EcoRV.cut(BamHI, EcoRI).pop(0)

pcaps = pCAPs.cut(PvuI).pop()

a = Assembly2([psu, pcaps], limit=200)

pYPK0 = a.circular_products[0]

pYPK0=pYPK0.synced(pCAPs)

pYPK0.name = "pYKP0"

