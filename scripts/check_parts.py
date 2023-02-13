import re
from pathlib import Path
from pathvalidate import ValidationError, validate_filename

text = Path("YeastPathwayKitSystematicNames.txt").read_text()

names = [ln for ln in text.splitlines() if "_" in ln]

datadirs = [
    "/home/bjorn/Desktop/mec@githb/YeastPathwayKit/sequences",
    "/home/bjorn/Desktop/mec@githb/YeastPathwayKitPrivate/sequences",]

vectors = []

for name in names:

    validate_filename(name, check_reserved=True)

    bb, *prts = name.split("_")

    vectors.append(f"{bb}")

    if "." in name:
        raise Exception("dot (.) not permitted in name")

    assert len(prts) % 2

    elms = [e.lower() for e in prts]

    if len(prts) != len(set(elms)):
        seen = set()
        dupl = []
        for i, e in enumerate(elms):
            if e.lower() in seen:
                dupl.append(elms.index(e.lower()))
                dupl.append(i)
            seen.add(e.lower())
        pattern = name
        for elm in set(prts[i] for i in dupl):
            pattern = pattern.replace(elm, len(elm)*"*")
        pattern = re.sub("[^*]", " ", pattern)
        raise ValueError(f"Duplicated elements: "
                         f"{' '.join(set(prts[i] for i in dupl))}\n"
                         f"{name}\n"
                         f"{pattern}\n")

    for i in range(0, len(prts)-2, 2):
        vectors.extend((f"pYPKa_Z_{prts[i]}",
                        f"pYPKa_A_{prts[i+1]}",
                        f"pYPKa_E_{prts[i+2]}"))

vectors = {k: None for k in vectors}.keys()

vectors_in_repo = []
for dd in datadirs:
    vectors_in_repo.extend(p.stem for p in Path(dd).rglob("*.gb"))

vectors_in_repo = {k: None for k in vectors_in_repo}.keys()

not_in_repo = set(vectors) - set(vectors_in_repo)

not_in_repo2 = set(v.lower() for v in vectors) - set(v.lower() for v in vectors_in_repo)

print(len(not_in_repo), len(not_in_repo2))

zzz = set([x[5:].lower() for x in vectors]) - set([v[5:].lower() for v in vectors_in_repo])


"""


pYPKa_E_GAS1
pYPKa_E_RPL4
pYPKa_E_RPS13


pYPKa_Z_RPL4
pYPKa_Z_RPS13
"""


from pydna.myprimers import PrimerList

p = PrimerList()

"""
pYPKa_A_AthMER
    pYPKa_A_AthfabD
    pYPKa_A_AthfabG
    pYPKa_A_Athfabz1
        pYPKa_A_Athfabz2
pYPKa_A_AthkasI
pYPKa_A_AthkasII
pYPKa_A_AthkasIII
pYPKa_A_BsfabHB
pYPKa_A_BsfabL
pYPKa_A_BslplJ
pYPKa_A_CaFATA
pYPKa_A_CaFATB
pYPKa_A_CvFATB1
pYPKa_A_EcAthacp1
pYPKa_A_EclplA
pYPKa_A_IAS1
pYPKa_A_IASA
pYPKa_A_LmfabH
pYPKa_A_McD9D
pYPKa_A_OAT
pYPKa_A_SafabH
pYPKa_A_SalplA
pYPKa_A_ScACC1
pYPKa_A_ScACC1S659AS1157A
pYPKa_A_SjACC1
pYPKa_A_UcFATB2
pYPKa_A_YlACC1
pYPKa_A_YlHMGL
pYPKa_A_YlIVD
pYPKa_A_YlMCCa
pYPKa_A_YlMCCb
pYPKa_A_YlMGH
pYPKa_A_YlbkdE1
pYPKa_A_YlbkdE1αΔMTS
pYPKa_A_YlbkdE1β
pYPKa_A_YlbkdE1βΔMTS
pYPKa_A_YlbkdE2
pYPKa_A_YlbkdE2ΔMTS
pYPKa_A_YlbkdE3
pYPKa_A_YlbkdE3ΔMTS




"""

"""

from pydna.sequence_picker import genbank_accession


x = genbank_accession(
    "GGATAAAAAAAAAAAGTTATTTTGAGTAGCTGATAAAGCGAGCTGGTGCCTATCATAGCCGGCTCAGACTTTTTATGAATTCACAGGCCAGCCCTGGCTATTCTTTTGCGTACTTTTAGTTCGATATATTTTCGCGGCTCGCGTTTTGTTTGCTTCTTATTTTACACTGAGTTTTCGTGCCGCAAACGTGGAGATGGGAAAAAGAAAAGTCGGGAAAATAATGAGAAATTTCTACTTTTGGTATTCCTCATACAGCCTGCGCGGTTTATTAGTAAAATACCCGATAATCCTCGAGGTTTGAAAAACTTTTCCCTCTACTACTGTTGACACGGATTTTTTTATTTAAGAGGAAAAGTCGTGGTTGTTTTCCTCGAACAAATTAGATATCCATAAATAGTTGTGTCGTTTTATTAAGCTATTTCAAAATCAGTTTTTATTTTTAAAGTCTGATAAAACAAAAACAACAAACACAGCTAAATCTCAACA")


y = genbank_accession(
"ATAATATTGAATATGCACTTTTACTATATTAATATCACGTCACACGACGCACAGTGAGAAGTGAAAAATTTTTTTTCAATCTGAAAAAAAAAAAAAAAAAAAAAAAAAATTTATATAAACGAATGGTATCTCCATCACATTTCTTTTAGCCTCGCAACTTGTACTTTTCATCACTTTTCTTGTAATTTAGCAATATCCCAAGAACAATCATCGAA")


z =  genbank_accession(
    
"TATAGTCTATCTTAAACACTAAACTACCTCCTATAATCATGTAGTGTACTTTAAACATTTTTTTATCTTCATAGCAATAATATAAGCCTTTTACCACCCATAAACCATAAAGTAGACCCAAACATTTTTAAAAAAATTTTACGTTATAATTTTTTCTTTGTCGTTTTTTCTGAGCGCGCAAAGTAGCGGTGAAATTTTGATACGAATGAGATTTCCACTTCTGTACAGATGGAAATTTATGTTGGCCGACATATATCACAGTCGTGATTGAATTAACAATTTCTTTCTCATTAATATTTATTCTAAACGGTTAACCACTAAATCAATCAACAACAATCAGTCAAA"
    )
"""
