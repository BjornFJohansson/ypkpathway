==============
pYPKa_E_PDC1tp
==============

Plan for the construction of E. coli vector `pYPKa_E_PDC1tp <./pYPKa_E_PDC1tp.txt>`_

Step 1 PCR of the insert
........................

PCR with primers pfw955 & prv955 and template `PDC1_template <./PDC1_template.txt>`_ results in 
a 968bp `PCR product <./PDC1.txt>`_


Primers annealing on template:
::

       5AGGGTAGCCTCCCCAT...ACAGTCAAATCAATCAAA3
                           |||||||||||||||||| tm 41.1 (dbd) 52.7
                          3TGTCAGTTTAGTTAGTTTAATTAAT5
 5TTAAATAGGGTAGCCTCCCCAT3
        |||||||||||||||| tm 49.1 (dbd) 61.3
       3TCCCATCGGAGGGGTA...TGTCAGTTTAGTTAGTTT5

Suggested PCR programs for Taq polymerase and for Polymerases with DNA binding domain:
::

 
 Taq (rate 30 nt/s) 35 cycles             |968bp
 95.0°C    |95.0°C                 |      |SantaLucia 1998
 |_________|_____          72.0°C  |72.0°C|SaltC 50mM
 | 03min00s|30s  \         ________|______|
 |         |      \ 50.0°C/ 0min29s| 5min |
 |         |       \_____/         |      |
 |         |         30s           |      |4-12°C

Step 2 Vector digestion and cloning
...................................

Clone the `PCR product <./PDC1.txt>`_ in `pYPKa <./pYPKa.txt>`_ digested 
with `EcoRV <http://rebase.neb.com/rebase/enz/EcoRV.html>`_ resulting in `pYPKa_E_PDC1tp <./pYPKa_E_PDC1tp.txt>`_


Step 3 Diagnostic PCR confirmation
..................................

Confirm the structure of the `pYPKa_E_PDC1tp <./pYPKa_E_PDC1tp.txt>`_ using primers 568, 342 and pfw955 
in a multiplex PCR reaction.

Expected PCR products sizes from 568, 342 and pfw955 (bp):

pYPKa with insert in correct orientation: 1684, 1653 |br|
pYPKa with insert in reverse orientation: 1684, 999 |br|
Empty pYPKa clone                       : 716 


.. |br| raw:: html

   <br />
