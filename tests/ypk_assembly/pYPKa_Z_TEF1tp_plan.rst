==============
pYPKa_Z_TEF1tp
==============

Plan for the construction of E. coli vector `pYPKa_Z_TEF1tp <./pYPKa_Z_TEF1tp.txt>`_

Step 1 PCR of the insert
........................

PCR with primers pfw579 & prv579 and template `TEF1_template <./TEF1_template.txt>`_ results in 
a 592bp `PCR product <./TEF1.txt>`_


Primers annealing on template:
::

       5ACAATGCATACTTTGTAC...TAATCTAAGTTTTAATTACAAA3
                             |||||||||||||||||||||| tm 39.1 (dbd) 49.3
                            3ATTAGATTCAAAATTAATGTTTaattaat5
 5ttaaatACAATGCATACTTTGTAC3
        |||||||||||||||||| tm 41.9 (dbd) 49.5
       3TGTTACGTATGAAACATG...ATTAGATTCAAAATTAATGTTT5

Suggested PCR programs for Taq polymerase and for Polymerases with DNA binding domain:
::

 
 Taq (rate 30 nt/s) 35 cycles             |592bp
 95.0°C    |95.0°C                 |      |SantaLucia 1998
 |_________|_____          72.0°C  |72.0°C|SaltC 50mM
 | 03min00s|30s  \         ________|______|
 |         |      \ 48.0°C/ 0min17s| 5min |
 |         |       \_____/         |      |
 |         |         30s           |      |4-12°C

Step 2 Vector digestion and cloning
...................................

Clone the `PCR product <./TEF1.txt>`_ in `pYPKa <./pYPKa.txt>`_ digested 
with `ZraI <http://rebase.neb.com/rebase/enz/ZraI.html>`_ resulting in `pYPKa_Z_TEF1tp <./pYPKa_Z_TEF1tp.txt>`_


Step 3 Diagnostic PCR confirmation
..................................

Confirm the structure of the `pYPKa_Z_TEF1tp <./pYPKa_Z_TEF1tp.txt>`_ using primers 577, 342 and pfw579 
in a multiplex PCR reaction.

Expected PCR products sizes from 577, 342 and pfw579 (bp):

pYPKa with insert in correct orientation: 1526, 1358 |br|
pYPKa with insert in reverse orientation: 1526, 760 |br|
Empty pYPKa clone                       : 934 


.. |br| raw:: html

   <br />
