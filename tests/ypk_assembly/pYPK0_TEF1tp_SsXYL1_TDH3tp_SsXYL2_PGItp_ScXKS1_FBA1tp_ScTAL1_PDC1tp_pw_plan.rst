======================================================================
pYPK0_TEF1tp_SsXYL1_TDH3tp_SsXYL2_PGItp_ScXKS1_FBA1tp_ScTAL1_PDC1tp_pw
======================================================================

Step 1 Prepare vector
.....................

Linearize `pYPKpw <./pYPKpw.txt>`_ with `EcoRV <http://rebase.neb.com/rebase/enz/EcoRV.html>`_
resulting in the `linearized vector <./pYPKpw_lin.txt>`_.

Step 2 tp-gene-tp PCR reactions
...............................

Perform the following 4 PCR reactions:

primers 577, 778 and `pYPK0_TEF1tp_SsXYL1_TDH3tp <./pYPK0_TEF1tp_SsXYL1_TDH3tp.txt>`__ => `2524bp_PCR_prod <./2524bp_PCR_prod.txt>`__ |br| |br|
primers 775, 778 and `pYPK0_TDH3tp_SsXYL2_PGItp <./pYPK0_TDH3tp_SsXYL2_PGItp.txt>`__ => `2924bp_PCR_prod <./2924bp_PCR_prod.txt>`__ |br| |br|
primers 775, 778 and `pYPK0_PGItp_ScXKS1_FBA1tp <./pYPK0_PGItp_ScXKS1_FBA1tp.txt>`__ => `3567bp_PCR_prod <./3567bp_PCR_prod.txt>`__ |br| |br|
primers 775, 578 and `pYPK0_FBA1tp_ScTAL1_PDC1tp <./pYPK0_FBA1tp_ScTAL1_PDC1tp.txt>`__ => `3009bp_PCR_prod <./3009bp_PCR_prod.txt>`__ |br| |br|



Step 3 Transformation and Assembly
..................................

Mix the DNA fragments and transform a S. cerevisiae ura3 mutant. The DNA fragments 
will be assembled by in-vivo homologous recombination:
::

  -|pYPKpw|124
 |         \/
 |         /\
 |         124|2524bp_PCR_prod|711
 |                             \/
 |                             /\
 |                             711|2924bp_PCR_prod|1013
 |                                                 \/
 |                                                 /\
 |                                                 1013|3567bp_PCR_prod|643
 |                                                                      \/
 |                                                                      /\
 |                                                                      643|3009bp_PCR_prod|242
 |                                                                                          \/
 |                                                                                          /\
 |                                                                                          242-
 |                                                                                             |
  ---------------------------------------------------------------------------------------------



.. |br| raw:: html

   <br />

