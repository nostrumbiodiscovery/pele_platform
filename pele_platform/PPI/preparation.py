import os
import subprocess

from pele_platform.constants import constants


def prepare_structure(protein_file, ligand_pdb, chain, remove_water=False, peleffy=False):
    
    to_remove = []

    # remove additional protein chains and all water molecules
    with open(protein_file, "r") as file:
        lines = file.readlines()

        for line in lines:
            if ((line.startswith("ATOM") or line.startswith("HETATM")) and line[21:22].strip() not in chain) or \
                    line.startswith("END") or line.startswith("TER") or line.startswith("CONECT"):
                to_remove.append(line)
            if remove_water:
                if (line.startswith("ATOM") or line.startswith("HETATM")) and line[17:20].strip() == "HOH":
                    to_remove.append(line)

        protein = [line for line in lines if line not in to_remove]
    
    # read in ligand file
    ligand = []
    with open(ligand_pdb, "r") as ligand_file:
        lines = ligand_file.readlines()
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                ligand.append(line)
    
    new_protein_file = os.path.basename(protein_file).replace(".pdb", "_prep.pdb")
    new_protein_file = os.path.abspath(new_protein_file)

    # join protein and ligand PDBs into new file
    with open(new_protein_file, "w+") as file:
        for line in protein:
            file.write(line)
        file.write("\n")
        for line in ligand:
            file.write(line)

    if peleffy:
        # Run through Schrodinger again to add CONECT lines.
        schrodinger_path = os.path.join(constants.SCHRODINGER, "utilities/prepwizard")
        command_pdb = f"{schrodinger_path} -nohtreat -noepik -noprotassign -noimpref -noccd -NOJOBID {protein_file} {new_protein_file}"
        subprocess.call(command_pdb.split(), shell=True)

    return new_protein_file
