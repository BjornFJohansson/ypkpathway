=============
pYPKa_Z_PGItp
=============

Plan for the construction of E. coli vector `pYPKa_Z_PGItp <./pYPKa_Z_PGItp.txt>`_

Step 1 PCR of the insert
........................

PCR with primers pfw1000 & prv1000 and template `PGI_template <./PGI_template.txt>`_ results in 
a 1013bp `PCR product <./PGI.txt>`_


Primers annealing on template:
::

       5AATTCAGTTTTCTGACTGA...CAAGATACCAGCCTAAAA3
                              |||||||||||||||||| tm 43.0 (dbd) 54.5
                             3GTTCTATGGTCGGATTTTAATTAAT5
 5TTAAATAATTCAGTTTTCTGACTGA3
        ||||||||||||||||||| tm 43.6 (dbd) 53.8
       3TTAAGTCAAAAGACTGACT...GTTCTATGGTCGGATTTT5

Suggested PCR programs for Taq polymerase and for Polymerases with DNA binding domain:
::

 
 Taq (rate 30 nt/s) 35 cycles             |1013bp
 95.0°C    |95.0°C                 |      |SantaLucia 1998
 |_________|_____          72.0°C  |72.0°C|SaltC 50mM
 | 03min00s|30s  \         ________|______|
 |         |      \ 51.0°C/ 0min30s| 5min |
 |         |       \_____/         |      |
 |         |         30s           |      |4-12°C

Step 2 Vector digestion and cloning
...................................

Clone the `PCR product <./PGI.txt>`_ in `pYPKa <./pYPKa.txt>`_ digested 
with `ZraI <http://rebase.neb.com/rebase/enz/ZraI.html>`_ resulting in `pYPKa_Z_PGItp <./pYPKa_Z_PGItp.txt>`_


Step 3 Diagnostic PCR confirmation
..................................

Confirm the structure of the `pYPKa_Z_PGItp <./pYPKa_Z_PGItp.txt>`_ using primers 577, 342 and pfw1000 
in a multiplex PCR reaction.

Expected PCR products sizes from 577, 342 and pfw1000 (bp):

pYPKa with insert in correct orientation: 1947, 1779 |br|
pYPKa with insert in reverse orientation: 1947, 1181 |br|
Empty pYPKa clone                       : 934 


.. |br| raw:: html

   <br />
