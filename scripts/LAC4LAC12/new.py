from ypkpathway import PathWay

workdir = "/home/bjorn/python_packages/ypkpathway/scripts/LAC4LAC12"

datafolders = (
    "/home/bjorn/Desktop/mec@githb/YeastPathwayKit/sequences",
    "/home/bjorn/Desktop/mec@githb/YeastPathwayKitPrivate/sequences",)

name = "pYPK0_PDC1tp_KlLAC4_PGI1tp_KlLAC12_TPI1tp.gb"

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

