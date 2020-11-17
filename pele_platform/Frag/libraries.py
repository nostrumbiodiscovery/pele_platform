import glob
import os
import subprocess

from pele_platform.constants import constants as cs

def growing_sites(fragment, user_bond):
    """
    Retrieves all possible growing sites (hydrogens) on the fragment. Takes PDB fragment file as input.
    Output - list of strings represeting sites, e.g. "benzene.pdb C6-H6 C1-H2"
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    bonds = []
    mol = Chem.MolFromPDBFile(fragment, removeHs=False)
    
    if mol:
        heavy_atoms = [a for a in mol.GetAtoms() if a.GetSymbol() != "H"]
        for a in heavy_atoms:
            hydrogens = [n for n in a.GetNeighbors() if n.GetSymbol() == "H"]
            at_name = a.GetMonomerInfo().GetName().strip()
            for h in hydrogens:
                h_name = h.GetMonomerInfo().GetName().strip()
                bonds.append("{} {} {}-{}".format(fragment, user_bond, at_name, h_name))

    return bonds


def extract_from_sdf(file_list, path):

    converted_mae = []
    output = []

    # convert all SDF to MAE
    schrodinger_path = os.path.join(cs.SCHRODINGER, "utilities/structconvert")
    command_mae = "{} -isd {} -omae {}"
    command_pdb = "{} -imae {} -opdb {}"

    for f in file_list:
        fout = os.path.splitext(os.path.basename(f))[0] + ".mae"
        fout_path = os.path.join(os.path.dirname(f), fout)
        try:
            command_mae = command_mae.format(schrodinger_path, f, fout_path)
            subprocess.call(command_mae.split())
            converted_mae.append(fout_path)
        except Exception as e:
            print(e)
    
    # convert all MAE to PDB
    for c in converted_mae:
        fout = c.replace(".mae",".pdb")
        try:
            command_pdb = command_pdb.format(schrodinger_path, c, fout)
            subprocess.call(command_pdb.split())
            os.remove(c)
        except Exception as e:
            print(e)
    
    pdb_pattern = os.path.splitext(converted_mae[0])
    converted_pdb = glob.glob(pdb_pattern[0]+"*"+".pdb")
    
    # ~~~ If it's stupid but it works (?), it isn't stupid. ~~~
    
    # read in PDB file created by Schrodinger, substitute residue name and add chain ID
    for c in converted_pdb:
        with open(c, "r") as fin:
            lines = fin.readlines()
            new_lines = []
            for line in lines:
                if line.startswith("HETATM") or line.startswith("ATOM"):
                    new_lines.append(line)
        
        new_lines = [l.replace("UNK", "GRW") for l in new_lines if "UNK" in l]
        new_lines = [l[:21]+"L"+l[22:] for l in new_lines]

        with open(c, "w") as fout:
            for line in new_lines:
                fout.write(line)
    return converted_pdb

def main(user_bond, frag_library):

    # get absolute path to the fragments library
    directory = os.path.dirname(os.path.abspath(__file__))
    path = frag_library if os.path.exists(frag_library) else os.path.join(directory, "Libraries", frag_library.strip())
    if not os.path.exists(path):
        raise OSError(f"File {frag_library} doesn't exist and is not one of our internal libraries. Please check the frag_library flag in input.yaml.")

    # get fragment files
    fragment_files = []
    extensions = ['*.pdb', '*.sdf']
    for e in extensions:
        fragment_files.extend(glob.glob(os.path.join(path, e.upper())))
        fragment_files.extend(glob.glob(os.path.join(path, e.lower())))

    # convert SDF to PDB, if necessary
    sdf_files = [elem for elem in fragment_files if ".sdf" in elem.lower()]
    all_files = [elem for elem in fragment_files if ".pdb" in elem.lower()]

    if sdf_files:
        all_files.extend(extract_from_sdf(sdf_files, path))

    # scan fragments for attachment points and write to input.conf
    bond_list = []
    for file in all_files:
        bond_list.extend(growing_sites(file, user_bond))

    output_name = "input.conf"
    with open(output_name, "w+") as conf_file:
        for line in bond_list:
            conf_file.write(line+"\n")

    return output_name
