from pathlib import Path
from pydna.readers import read

def process_pathway_name(name, TU_bb_name="pYPK0"):
    
    backbone, *items = name.split("_")
    
    assert len(items) % 2 == 1 # Uneven number of elements
    
    assert len(items) >= 3 # At least one TU
    
    backbone_vector_name = Path(f"{backbone}.gb")
    
    tu_vector_names = []
    
    for i in range(0,len(items)-1, 2):
        p, g, t = items[i:i+3]
        tu_vector_names.append(Path(f"{TU_bb_name}_{p}_{g}_{t}.gb"))
        
    return backbone_vector_name, tu_vector_names
  

# name = "pTA1_TDH3_ScATF1_PGI1"
name = "pTA1_TDH3_ScATF1_PGI1_CTT1_TEF1"

backbone_vector_name, tu_vector_names = process_pathway_name(name)

def assemble_tu_vectors(tu_vector_names, userpath=Path("GeneticAssembly")):
    for tu_name in tu_vector_names:
        
        




userpath = Path("GeneticAssembly") 



promoters, genes, terminators = items[0:-1:2], items[1::2], items[2::2]



promoter_vector_names = [Path(f"pYPKa_Z_{p}.gb") for p in promoters]

gene_vector_names = [Path(f"pYPKa_A_{p}.gb") for p in genes]

terminator_vector_names = [Path(f"pYPKa_E_{p}.gb") for p in terminators]

tu_vector_names = []

for p, g, t in zip(promoters, genes, terminators):
    tu_vector_names.append(Path(f"pYPK0_{p}_{g}_{t}.gb"))

# userpath = Path("/content/gdrive/My Drive/GeneticAssembly")


userpath.mkdir(parents=True, exist_ok=True)

ypkpath = Path("YeastPathwayKit/sequences")

pYPKa_As = []

for pth in ypkpath.glob("pYPKa_Z_*.gb"):
    pYPKa_Zs.append(read(pth))
    
for pth in userpath.glob("pYPKa_A_*.gb"):
    pYPKa_Zs.append(read(pth))
    
for pth in ypkpath.glob("pYPKa_E_*.gb"):
    pYPKa_Es.append(read(pth))

