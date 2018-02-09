import os
from schrodinger import structure as st
import StringIO
import MSM_PELE.Helpers.helpers as hp

def build_complexes(ligands, receptor):
    """
       Desciption: From each structure retrieve
       a .mae file, the complex with the receptor
       and a the ligand's residue name.

       Output:
            struct_files_mae: list of ligand.mae files
            complex: list of receptor+ligand pdb
            residues: list of residues of ligands
    """

    
    for structure in st.StructureReader(ligands):
        for residue in structure.residue:
            res = residue.pdbres.strip()
        str_name = "{}".format(res)
        try:
            structure.write(str_name + ".mae")
        except ValueError:
            str_name = "{}".format(res)
        finally:
            structure.write(str_name + ".mae")
            structure_mae = "{}.mae".format(str_name)
    return structure_mae, res


def retrieve_receptor(system, residue):
    """
    This function returns receptor of the complex of interest.

    :param complex: system format pdb

    :output: receptor text
    """
    ligand = os.path.abspath("lig.pdb")
    with open(system, 'r') as pdb_file:
        receptor_text = [line for line in pdb_file if line.startswith("ATOM")]
    with open(system, 'r') as pdb_file:
        ligand_text = [line for line in pdb_file if line[17:20].strip() == residue]
    if not receptor_text  or not ligand_text:
        raise ValueError("Something went wrong when extracting the ligand. Ligand must be a HETATOM")
    with open(ligand, "w") as fout:
	fout.write("".join(ligand_text))

    return "".join(receptor_text), ligand 
