from ypkpathway import PathWay





pws = [c2 for (c1,c2) in [line.split() for line in 
"""\
pEc1    pTA1_PDC1_EcfabH_TEF1_EcfabD_FBA1_EcfabG_RPL22A_EcacpP_TDH3_EcfabF_UTR2_EcfabB_TPI1_EcfabA_PMP3_EcfabZ_ENO2_EcfabI_RPL5_EctesA_RPL16A_EcacpS_RPL17A
pEc2    pTA1_PDC1_EcfabH_TEF1_EcfabD_FBA1_EcfabG_RPL22A_EcacpP_TDH3_EcfabF_UTR2_EcfabB_TPI1_EcfabA_PMP3_EcfabZ_ENO2_EcfabI_RPL5_EctesB_RPL16A_EcacpS_RPL17A
pEc3    pTA1_PDC1_EcfabH_TEF1_EcfabD_FBA1_EcfabG_RPL22A_EcacpP_TDH3_EcfabF_UTR2_EcfabB_TPI1_EcfabA_PMP3_EcfabZ_ENO2_EcfabI_RPL5_EctesA_RPL16A_EctesB_RPL17A_EcacpS_TMA19
pEcBs   pTA1_PDC1_EcfabH_TEF1_EcfabD_FBA1_EcfabG_RPL22A_EcacpP_TDH3_EcfabF_UTR2_EcfabB_TPI1_EcfabA_PMP3_EcfabZ_ENO2_BsfabL_RPL5_EctesA_RPL16A_EctesB_RPL17A_EcacpS_TMA19
pAt19   pTA1_PDC1_AthkasIII_TEF1_AthfabD_FBA1_AthfabG_RPL22A_EcAthacp1_TDH3_AthkasII_UTR2_AthkasI_TPI1_Athfabz1_PMP3_Athfabz2_ENO2_Athmod1_RPL5_AthfatA1_RPL16A_EcacpS_RPL17A
pAt20   pTA1_PDC1_AthkasIII_TEF1_AthfabD_FBA1_AthfabG_RPL22A_EcAthacp1_TDH3_AthkasII_UTR2_AthkasI_TPI1_Athfabz1_PMP3_Athfabz2_ENO2_Athmod1_RPL5_AthfatB_RPL16A_EcacpS_RPL17A
pAt21   pTA1_PDC1_AthkasIII_TEF1_AthfabD_FBA1_AthfabG_RPL22A_EcAthacp1_TDH3_AthkasII_UTR2_AthkasI_TPI1_Athfabz1_PMP3_Athfabz2_ENO2_Athmod1_RPL5_AthfatA1_RPL16A_AthfatB_RPL17A_EcacpS_TMA19
pAtMER  pTA1_PDC1_AthkasIII_TEF1_AthfabD_FBA1_AthfabG_RPL22A_EcAthacp1_TDH3_AthkasII_UTR2_AthkasI_TPI1_Athfabz1_PMP3_Athfabz2_ENO2_AthMER_RPL5_AthfatA1_RPL16A_AthfatB_RPL17A_EcacpS_TMA19
""".splitlines()]]


name = pws[3]

workdir = "/home/bjorn/Desktop/python_packages/ypkpathway/scripts/humberto"

datafolders = (
    "/home/bjorn/Desktop/mec@githb/YeastPathwayKit/sequences",
    "/home/bjorn/Desktop/mec@githb/YeastPathwayKitPrivate/sequences",)
    
pw = PathWay(name, workdir, datafolders, tu_backbone="pYPK0")

pw.find_sequence_files()
pw.copy_sequence_files()
pw.read_sequence_files()

pw.find_primers()
pw.copy_primers()
pw.read_primers()
pw.find_enzymes()

pw.execute_transcriptional_unit_notebooks()
pw.find_sequence_files()
pw.read_sequence_files()
pw.pcr()
pw.execute_nb()

