from pathlib import Path
from pydna.readers import read
import importlib
import jupytext
from Bio.Restriction import ZraI, AjiI, EcoRV


def cloned(vector, enzyme, candidate):
    if len(candidate) <= len(vector):
        return 0
    candidate2 = str(candidate.seq).lower()*2
    for linear_vector in vector.cut(enzyme):
        if str(linear_vector.seq).lower() in candidate2:
            return len(candidate) - len(vector)
    return 0


def df(name):
    f = importlib.resources.files("ypkpathway")/Path("data")/name
    return f.read_text()


def load_files(pth="GeneticAssembly"):
    pth = Path(pth)
    pth.mkdir(parents=True, exist_ok=True)
    return list(pth.glob("*.gb"))


def process_pathway_name(name, TU_bb_name="pYPK0"):
    
    backbone, *items = name.split("_")
    
    assert len(items) % 2 == 1 # Uneven number of elements
    
    assert len(items) >= 3 # At least one TU
    
    backbone_vector_name = Path(f"{backbone}.gb")

    tu_vector_names = []

    for i in range(0, len(items)-1, 2):
        p, g, t = items[i:i+3]
        tu_vector_names.append(Path(f"{TU_bb_name}_{p}_{g}_{t}.gb"))

    return backbone_vector_name, tu_vector_names


def assemble_tu_vectors(tu_vector_names,
                        usergbfiles,
                        kitfgbfiles):
    files = {}
    for tu_name in tu_vector_names:
        if tu_name.name in [f.name for f in usergbfiles]:
            files[tu_name] = tu_name.read_text()
        else:
            files[tu_name] = TU_assembly(tu_name)
    return files
    

def TU_assembly(tu_name):
    
    tu_name = Path(tu_name)

    files = {}

    bb_name, promoter_name, gene_name, terminator_name = tu_name.stem.split("_")

    promoter_vector_path = kit_folder/Path(f"pYPKa_Z_{promoter_name}.gb")
    gene_vector_path = user_folder/Path(f"pYPKa_A_{gene_name}.gb")
    terminator_vector_path = kit_folder/Path(f"pYPKa_E_{terminator_name}.gb")

    mdnbtemplate = df("nb_template_pYPK0_tp_gene_tp.md")
    
    formattedmd = mdnbtemplate.format(tpz=promoter_name,
                                      gene=gene_name,
                                      tpe=terminator_name)
    
    ipynb = jupytext.reads(formattedmd, fmt='md')
    
    files[tu_name.with_suffix(".md")] = jupytext.writes(ipynb, fmt='ipynb')
    
    return files


def execute_notebooks(files, log, working_dir="ypkassembly", print=print):

    try:
        pathlib.Path(working_dir).mkdir(parents=True)
    except OSError as exception:
        if exception.errno == errno.EEXIST:
            print("The {} directory already exists! Please delete or choose another name.".format(working_dir))
        else:
            print("The {} directory could not be created".format(working_dir))
        return None, None


    cwd = os.getcwd()

    msg = "created subdirectory {}\n".format(working_dir)

    print(msg)

    log+=msg

    os.chdir(working_dir)

    msg = "saving files sequence files and images.\n"
    print(msg)
    log+=msg

    notebooks = []
    obj = notedown.MarkdownReader()

    for name, content in files.items():
        if name.endswith(".md"):
            newname = os.path.splitext(name)[0]+".ipynb"
            msg = "\nsaving: "+newname
            print(msg)
            log+=msg
            nb = nbformat.write(obj.to_notebook(content), newname)
            notebooks.append(newname)
        else:
            mode = {True:"wb", False:"w"}[hasattr(content, "decode")]
            with open(name, mode) as f:
                f.write(content)

    del files

    ep = ExecutePreprocessor()
    ep.timeout = 1200              # seconds
    ep.interrupt_on_timeout = True

    print("\n")
    log+="\n"

    msg = "\nexecuting pYPKa notebooks..\n"
    print(msg)
    log+=msg

    shell = InteractiveShell.instance()

    g={}
    l={}
    resources = {}
    
    for name in notebooks[:-1]:
        msg = "\nexecuting: "+name
        print(msg)
        log+=msg
        with open(name, 'r', encoding='utf-8') as f:
            nb = nbformat.read(f, 4)
        nb_executed, resources = ep.preprocess(nb,
                                               resources=resources)
        for cell in nb.cells:
            if cell.cell_type == 'code':
                code = shell.input_transformer_manager.transform_cell(cell.source)
                exec(code, g, l)
        nbformat.write(nb, name)
        g={}
        l={}

    os.chdir(cwd)

    fl = FileLink(os.path.join(working_dir, "pw.ipynb"))

    return fl, log


if __name__ == "__main__":

    # name = "pTA1_TDH3_ScATF1_PGI1"
    
    name = "pTA1_TDH3_ScATF1_PGI1_ScCTT1_TEF1"
    
    name = "pTA1_PDC1_EcfabH_TEF1_EcfabD_FBA1_EcfabG_RPL22A_EcacpP_TDH3_EcfabF_UTR2_EcfabB_TPI1_EcfabA_PMP3_EcfabZ_ENO2_Athmod1_RPL5_AthfatA1_RPL16A_AthfatB_RPL17A_EcacpH_RPL16B_EcacpS_TMA19"
    
    backbone_vector_name, tu_vector_names = process_pathway_name(name)
    
    user_folder = Path("GeneticAssembly")
    kit_folder = Path("YeastPathwayKit/sequences")
    usergbfiles = load_files(user_folder)
    kitfgbfiles = load_files(kit_folder)
    
    x = assemble_tu_vectors(tu_vector_names, usergbfiles, kitfgbfiles)
    
    
