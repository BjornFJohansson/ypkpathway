==========================
pYPK0_FBA1tp_ScTAL1_PDC1tp
==========================

Step 1 Prepare vector
.....................

Linearize `pYPKpw <./pYPKpw.txt>`_ with `EcoRV <http://rebase.neb.com/rebase/enz/EcoRV.html>`_
resulting in the `linearized vector <./pYPKpw_lin.txt>`_.

Step 2 PCR of first tp
......................

Carry out a PCR with primers 577, 567 and template `pYPKa_Z_FBA1tp <./pYPKa_Z_FBA1tp.txt>`_ resulting in 
the PCR product `861bp_PCR_prod <./pYPKa_Z_FBA1tp_pcr_prd.txt>`_      |br|   
::

 5GTTCTGATCCTCGAGCATCTTAAGAATTC...CTCACTAGTGACCTGCAGCCGAC3
                                  ||||||||||||||||||||||| tm 59.0 (dbd) 70.7
                                 3gagtgatcactggacgtcggcTG5
 5gttctgatcctcgagcatcttaagaattc3
  ||||||||||||||||||||||||||||| tm 56.1 (dbd) 69.4
 3CAAGACTAGGAGCTCGTAGAATTCTTAAG...GAGTGATCACTGGACGTCGGCTG5

 
 Taq (rate 30 nt/s) 35 cycles             |861bp
 95.0°C    |95.0°C                 |      |SantaLucia 1998
 |_________|_____          72.0°C  |72.0°C|SaltC 50mM
 | 03min00s|30s  \         ________|______|
 |         |      \ 56.0°C/ 0min25s| 5min |
 |         |       \_____/         |      |
 |         |         30s           |      |4-12°C

Step 3 Gene PCR
...............

Carry out a PCR with primers 468, 467 and template `pYPKa_A_ScTAL1 <./pYPKa_A_ScTAL1.txt>`_ resulting in 
the PCR product `1097bp_PCR_prod <./pYPKa_A_ScTAL1_pcr_prd.txt>`_     |br|   
::

 5GTCGAGGAACGCCAGGTTGCCCACT...TCTGTGCAGACAAACGCATCAGGAT3
                              ||||||||||||||||||||||||| tm 59.0 (dbd) 73.8
                             3agacacgtctgtttgcgtagtcctaAATTTA5
 5gtcgaggaacgccaggttgcccact3
  ||||||||||||||||||||||||| tm 64.8 (dbd) 79.7
 3CAGCTCCTTGCGGTCCAACGGGTGA...AGACACGTCTGTTTGCGTAGTCCTA5

 
 Taq (rate 30 nt/s) 35 cycles             |1097bp
 95.0°C    |95.0°C                 |      |SantaLucia 1998
 |_________|_____          72.0°C  |72.0°C|SaltC 50mM
 | 03min00s|30s  \         ________|______|
 |         |      \ 58.0°C/ 0min32s| 5min |
 |         |       \_____/         |      |
 |         |         30s           |      |4-12°C

Step 4 PCR of last tp
.....................

Carry out a PCR with primers 568, 578 and template `pYPKa_E_PDC1tp <./pYPKa_E_PDC1tp.txt>`_ resulting in 
the PCR product `1294bp_PCR_prod <./pYPKa_E_PDC1tp_pcr_prd.txt>`_      |br|   
::

 5GTGCCATCTGTGCAGACAAACG...ACTTATGAATGTGGCAATGAGACAAGAAC3
                           ||||||||||||||||||||||||||||| tm 56.5 (dbd) 69.5
                          3tgaatacttacaccgttactctgttcttg5
 5GTGCcatctgtgcagacaaacg3
  |||||||||||||||||||||| tm 57.1 (dbd) 71.5
 3CACGGTAGACACGTCTGTTTGC...TGAATACTTACACCGTTACTCTGTTCTTG5

 
 Taq (rate 30 nt/s) 35 cycles             |1294bp
 95.0°C    |95.0°C                 |      |SantaLucia 1998
 |_________|_____          72.0°C  |72.0°C|SaltC 50mM
 | 03min00s|30s  \         ________|______|
 |         |      \ 55.0°C/ 0min38s| 5min |
 |         |       \_____/         |      |
 |         |         30s           |      |4-12°C


Step 5 Yeast transformation
...........................

Mix the four linear DNA fragments and transform a Saccharomyces cerevisiae ura3 mutant with the mixture.
The fragments will be assembled by in-vivo homologous recombination:

::

  -|pYPKpw|124
 |         \/
 |         /\
 |         124|861bp_PCR_prod|50
 |                            \/
 |                            /\
 |                            50|1097bp_PCR_prod|37
 |                                               \/
 |                                               /\
 |                                               37|1294bp_PCR_prod|242
 |                                                                  \/
 |                                                                  /\
 |                                                                  242-
 |                                                                     |
  ---------------------------------------------------------------------



Step 6 Diagnostic PCR confirmation
..................................

First tp and gene
+++++++++++++++++

PCR using primers 577 & 467 |br|     

PCR products (bp)

    Correct          : 1908 |br|
    Missing first tp : 1265 |br|
    Missing gene     : 898 |br|
    Missing both     : 255 |br|

Gene and last tp
++++++++++++++++

PCR using primers 468 & 578 |br| 

PCR products (bp)

    Correct         : 2354 |br|
    Missing gene    : 1344 |br|
    Missing last tp : 1386 |br|
    Missing both    : 376 |br|

.. |br| raw:: html

   <br />

