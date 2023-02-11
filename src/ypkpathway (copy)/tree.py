#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# https://anytree.readthedocs.io/en/latest/index.html

from anytree.render import ContRoundStyle

from anytree import Node, RenderTree
root = Node("pw")
s0 = Node("sub0", parent=root)
s0b = Node("sub0B", parent=s0)
s0a = Node("sub0A", parent=s0)
s1 = Node("sub1", parent=root)

print(RenderTree(root, style=ContRoundStyle()))

# pYPK0_PDC1tp_KlLAC4_PGI1tp_KlLAC12_TPI1tp.gb



root

print("""
pYPK0_PDC1_EcfabH_TEF1.gb 🔴
pYPKa_Z_PDC1	🟢✅😃
pYPKa_A_EcfabH	🟢✅😡
pYPKa_E_TEF1	🟢🚫
""")






