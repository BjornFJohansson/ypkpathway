==============
pYPKa_A_SsXYL1
==============

Plan for the construction of E. coli vector `pYPKa_A_SsXYL1 <./pYPKa_A_SsXYL1.txt>`_

Step 1 PCR of the insert
........................

PCR with primers pfw957 & prv957 and template `SsXYL1_template <./SsXYL1_template.txt>`_ results in 
a 959bp `PCR product <./SsXYL1.txt>`_


Primers annealing on template:
::

   5ATGCCTTCTATTAAGTTGAA...AAGATTCCTATCTTCGTCTAA3
                           ||||||||||||||||||||| tm 45.3 (dbd) 55.5
                          3TTCTAAGGATAGAAGCAGATT5
 5aaATGCCTTCTATTAAGTTGAA3
    |||||||||||||||||||| tm 43.8 (dbd) 54.9
   3TACGGAAGATAATTCAACTT...TTCTAAGGATAGAAGCAGATT5

Suggested PCR programs for Taq polymerase and for Polymerases with DNA binding domain:
::

 
 Taq (rate 30 nt/s) 35 cycles             |959bp
 95.0°C    |95.0°C                 |      |SantaLucia 1998
 |_________|_____          72.0°C  |72.0°C|SaltC 50mM
 | 03min00s|30s  \         ________|______|
 |         |      \ 54.0°C/ 0min28s| 5min |
 |         |       \_____/         |      |
 |         |         30s           |      |4-12°C

Step 2 Vector digestion and cloning
...................................

Clone the `PCR product <./SsXYL1.txt>`_ in `pYPKa <./pYPKa.txt>`_ digested 
with `AjiI <http://rebase.neb.com/rebase/enz/AjiI.html>`_ resulting in `pYPKa_A_SsXYL1 <./pYPKa_A_SsXYL1.txt>`_


Step 3 Diagnostic PCR confirmation
..................................

Confirm the structure of the `pYPKa_A_SsXYL1 <./pYPKa_A_SsXYL1.txt>`_ using primers 468, 342 and pfw957 
in a multiplex PCR reaction.

Expected PCR products sizes from 468, 342 and pfw957 (bp):

pYPKa with insert in correct orientation: 1725, 1675 |br|
pYPKa with insert in reverse orientation: 1725, 1009 |br|
Empty pYPKa clone                       : 766 


.. |br| raw:: html

   <br />
