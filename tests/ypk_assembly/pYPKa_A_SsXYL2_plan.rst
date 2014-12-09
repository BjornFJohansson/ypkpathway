==============
pYPKa_A_SsXYL2
==============

Plan for the construction of E. coli vector `pYPKa_A_SsXYL2 <./pYPKa_A_SsXYL2.txt>`_

Step 1 PCR of the insert
........................

PCR with primers pfw1092 & prv1092 and template `SsXYL2_template <./SsXYL2_template.txt>`_ results in 
a 1094bp `PCR product <./SsXYL2.txt>`_


Primers annealing on template:
::

   5ATGACTGCTAACCCTTC...TGACGGCCCTGAGTAA3
                        |||||||||||||||| tm 48.0 (dbd) 60.8
                       3ACTGCCGGGACTCATT5
 5aaATGACTGCTAACCCTTC3
    ||||||||||||||||| tm 44.5 (dbd) 54.0
   3TACTGACGATTGGGAAG...ACTGCCGGGACTCATT5

Suggested PCR programs for Taq polymerase and for Polymerases with DNA binding domain:
::

 
 Taq (rate 30 nt/s) 35 cycles             |1094bp
 95.0°C    |95.0°C                 |      |SantaLucia 1998
 |_________|_____          72.0°C  |72.0°C|SaltC 50mM
 | 03min00s|30s  \         ________|______|
 |         |      \ 55.0°C/ 0min32s| 5min |
 |         |       \_____/         |      |
 |         |         30s           |      |4-12°C

Step 2 Vector digestion and cloning
...................................

Clone the `PCR product <./SsXYL2.txt>`_ in `pYPKa <./pYPKa.txt>`_ digested 
with `AjiI <http://rebase.neb.com/rebase/enz/AjiI.html>`_ resulting in `pYPKa_A_SsXYL2 <./pYPKa_A_SsXYL2.txt>`_


Step 3 Diagnostic PCR confirmation
..................................

Confirm the structure of the `pYPKa_A_SsXYL2 <./pYPKa_A_SsXYL2.txt>`_ using primers 468, 342 and pfw1092 
in a multiplex PCR reaction.

Expected PCR products sizes from 468, 342 and pfw1092 (bp):

pYPKa with insert in correct orientation: 1860, 1810 |br|
pYPKa with insert in reverse orientation: 1860, 1144 |br|
Empty pYPKa clone                       : 766 


.. |br| raw:: html

   <br />
