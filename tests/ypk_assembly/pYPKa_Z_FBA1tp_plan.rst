==============
pYPKa_Z_FBA1tp
==============

Plan for the construction of E. coli vector `pYPKa_Z_FBA1tp <./pYPKa_Z_FBA1tp.txt>`_

Step 1 PCR of the insert
........................

PCR with primers pfw630 & prv630 and template `FBA1_template <./FBA1_template.txt>`_ results in 
a 643bp `PCR product <./FBA1.txt>`_


Primers annealing on template:
::

       5ATAACAATACTGACAGTACTAAA...ACCAAGTAATACATATTCAAA3
                                  ||||||||||||||||||||| tm 42.2 (dbd) 52.2
                                 3TGGTTCATTATGTATAAGTTTaattaat5
 5ttaaatATAACAATACTGACAGTACTAAA3
        ||||||||||||||||||||||| tm 44.9 (dbd) 51.0
       3TATTGTTATGACTGTCATGATTT...TGGTTCATTATGTATAAGTTT5

Suggested PCR programs for Taq polymerase and for Polymerases with DNA binding domain:
::

 
 Taq (rate 30 nt/s) 35 cycles             |643bp
 95.0°C    |95.0°C                 |      |SantaLucia 1998
 |_________|_____          72.0°C  |72.0°C|SaltC 50mM
 | 03min00s|30s  \         ________|______|
 |         |      \ 49.0°C/ 0min19s| 5min |
 |         |       \_____/         |      |
 |         |         30s           |      |4-12°C

Step 2 Vector digestion and cloning
...................................

Clone the `PCR product <./FBA1.txt>`_ in `pYPKa <./pYPKa.txt>`_ digested 
with `ZraI <http://rebase.neb.com/rebase/enz/ZraI.html>`_ resulting in `pYPKa_Z_FBA1tp <./pYPKa_Z_FBA1tp.txt>`_


Step 3 Diagnostic PCR confirmation
..................................

Confirm the structure of the `pYPKa_Z_FBA1tp <./pYPKa_Z_FBA1tp.txt>`_ using primers 577, 342 and pfw630 
in a multiplex PCR reaction.

Expected PCR products sizes from 577, 342 and pfw630 (bp):

pYPKa with insert in correct orientation: 1577, 1409 |br|
pYPKa with insert in reverse orientation: 1577, 811 |br|
Empty pYPKa clone                       : 934 


.. |br| raw:: html

   <br />
