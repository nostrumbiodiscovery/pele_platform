import os
from schrodinger import structure as st
import StringIO

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

    residues = []
    structures_mae = []
    structures_pdb = []

    for structure in st.StructureReader(ligands):
	    residues.append([residue for residue in structure.residue][0].pdbres.strip())

    for i, (structure, res) in enumerate(zip(st.StructureReader(ligands), residues)):

        str_name = "{}".format(res)
        if str_name in '\t'.join(structures_mae):
            str_name = "ligand_{:d}".format(i)
        try:
            structure.write(str_name + ".mae")
            structure.write(str_name + ".pdb")
        except ValueError:
            str_name = "ligand_{:d}".format(i)
        finally:
            structure.write(str_name + ".mae")
            structure.write(str_name + ".pdb")
            structures_mae.append("{}.mae".format(str_name))
            structures_pdb.append("{}.pdb".format(str_name))

    struct_files_mae = [os.path.abspath(structure) for structure in structures_mae]
    struct_files_pdb = [os.path.abspath(structure) for structure in structures_pdb]

    complexes = systems(receptor, struct_files_pdb)

    return complexes, struct_files_mae, residues


def systems(receptor, structs):
    """
       Description: From one ligand structure build the 
       complex with its receptor.

       Output:
            self.struct_files_pdb: List of complexes.
    """
    for struct in structs:
        with open(struct, "r") as pdb:
            clean_lines = [line for line in pdb.readlines() if line.startswith("HETATM")]
        system = receptor + "TER\n" + "".join(clean_lines) + "END\n"
        with open(struct, "w") as pdb_out:
            pdb_out.write(system)
    return structs


def retrieve_receptor(system, residue):
    """
    This function returns receptor of the complex of interest.

    :param complex: system format pdb

    :output: receptor text
    """
     
    ligand = StringIO.StringIO()
    
    with open(system, 'r') as pdb_file:
        receptor_text = [line for line in pdb_file if line.startswith("ATOM")]
        ligand_text = [line for line in pdb_file if line[17:20].strip() == residue]

    if receptor_text == "" or ligand_text == "" :
        raise ValueError("Something went wrong when extracting the ligand. Ligand must be a HETATOM")
   
    ligand.write("".join(ligand_text))
    
    return "".join(receptor_text), ligand
