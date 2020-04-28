import os


def prepare_structure(protein_file, ligand_pdb, chain):
    
    to_remove = []
    
    with open(protein_file, "r") as file:
        lines = file.readlines()
        for line in lines:
            if ((line.startswith("ATOM") or line.startswith("HETATM")) and line[21:22].strip() not in chain) or line.startswith("END"):
                to_remove.append(line)
        protein = [line for line in lines if line not in to_remove]
    
    # read in ligand file
    ligand = []
    print(ligand_pdb)
    with open(ligand_pdb, "r") as ligand_file:
        lines = ligand_file.readlines()
        for line in lines:
            print(line)
            if line.startswith("ATOM") or line.startswith("HETATM"):
                ligand.append(line)
    
    new_protein_file = os.path.basename(protein_file).replace(".pdb", "_prep.pdb")
    new_protein_file = os.path.abspath(new_protein_file)
    print(new_protein_file)
    
    # join protein and ligand pdbs into new file
    with open(new_protein_file, "w+") as file:
        for line in protein:
            file.write(line)
        file.write("\n")
        for line in ligand:
            file.write(line)
    
    return new_protein_file
