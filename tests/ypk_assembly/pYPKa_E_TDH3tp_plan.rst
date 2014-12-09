==============
pYPKa_E_TDH3tp
==============

Plan for the construction of E. coli vector `pYPKa_E_TDH3tp <./pYPKa_E_TDH3tp.txt>`_

Step 1 PCR of the insert
........................

PCR with primers pfw698 & prv698 and template `TDH3_template <./TDH3_template.txt>`_ results in 
a 711bp `PCR product <./TDH3.txt>`_


Primers annealing on template:
::

       5ATAAAAAACACGCTTTTTC...AACACACATAAACAAACAAA3
                              |||||||||||||||||||| tm 44.1 (dbd) 54.7
                             3TTGTGTGTATTTGTTTGTTTaattaat5
 5ttaaatATAAAAAACACGCTTTTTC3
        ||||||||||||||||||| tm 42.7 (dbd) 55.3
       3TATTTTTTGTGCGAAAAAG...TTGTGTGTATTTGTTTGTTT5

Suggested PCR programs for Taq polymerase and for Polymerases with DNA binding domain:
::

 
 Taq (rate 30 nt/s) 35 cycles             |711bp
 95.0°C    |95.0°C                 |      |SantaLucia 1998
 |_________|_____          72.0°C  |72.0°C|SaltC 50mM
 | 03min00s|30s  \         ________|______|
 |         |      \ 50.0°C/ 0min21s| 5min |
 |         |       \_____/         |      |
 |         |         30s           |      |4-12°C

Step 2 Vector digestion and cloning
...................................

Clone the `PCR product <./TDH3.txt>`_ in `pYPKa <./pYPKa.txt>`_ digested 
with `EcoRV <http://rebase.neb.com/rebase/enz/EcoRV.html>`_ resulting in `pYPKa_E_TDH3tp <./pYPKa_E_TDH3tp.txt>`_


Step 3 Diagnostic PCR confirmation
..................................

Confirm the structure of the `pYPKa_E_TDH3tp <./pYPKa_E_TDH3tp.txt>`_ using primers 568, 342 and pfw698 
in a multiplex PCR reaction.

Expected PCR products sizes from 568, 342 and pfw698 (bp):

pYPKa with insert in correct orientation: 1427, 1396 |br|
pYPKa with insert in reverse orientation: 1427, 742 |br|
Empty pYPKa clone                       : 716 


.. |br| raw:: html

   <br />
