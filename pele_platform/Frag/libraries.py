from rdkit import Chem
import yaml, glob, argparse, os
from string import Template

InBond = "$BOND"


def growing_sites(fragment):
    """
    Function that obtains the possible starting points to grow a ligand. 
    Entrance: File to obtain the possible bonds of a protein to commence the growth
    Exit: List of strings which represent possible starting bonds
    """
    bonds = []
        
    mol = Chem.MolFromPDBFile(fragment, removeHs=False)
    if mol:
        heavy_atoms = [a for a in mol.GetAtoms() if a.GetSymbol() != "H"]
        for a in heavy_atoms:
            hydrogens = [n for n in a.GetNeighbors() if n.GetSymbol() == "H"]
            at_name = a.GetMonomerInfo().GetName().strip()
            for h in hydrogens:
                h_name = h.GetMonomerInfo().GetName().strip()
                bonds.append("{}-{}".format(at_name, h_name)) 
    return bonds


def input_file(fragment_pdb, core_atom, growing_site, output="input.conf", createFile=True):
    """
    Function that creates the input.conf file to for the FragPELE module.
    Arguments:
        -Fragment_pdb: The file of the protein to be grown.
        -Core_atom: The starting bond of the fragment to grown.
        -Growing_site: list of bonds where the fragment growing will happen.
        -Output: file name to save the configuration, default input.conf.
        -CreateFile: Parameter used to create a new input.conf file or append to an existing one. 
    """
    if createFile:
        with open(output, "w+") as file:
            output = write_File(file, False, fragment_pdb, core_atom, growing_site, output)
    else:
        with open(output, "a+") as file:
            output = write_File(file, True, fragment_pdb, core_atom, growing_site, output)
    return output


def write_File(file, isAppending, fragment_pdb, core_atom, growing_site, output):
    StrForFile = ""
    for bond in growing_site:
        if isAppending:
            StrForFile += ("\n" + fragment_pdb + " $BOND " + bond + "\n")
            isAppending = False
        else:
            StrForFile += (fragment_pdb + " $BOND " + bond + "\n")
    s = Template(StrForFile)
    s = s.safe_substitute(BOND=core_atom)
    file.write(s[:-1])
    return file

def main(user_bond, frag_library):
    directory = os.path.dirname(__file__)
    path = frag_library if os.path.exists(frag_library) else os.path.join(directory, "Libraries", frag_library.strip(), "Core,*.pdb")
    
    if not os.path.exists(path):
        raise OSError("Cannot find the library.")
    
    path = os.path.abspath(path)
    firstTime = True
    fragment_files = glob.glob(os.path.join(path, "*.pdb"))
    print("Retrieved {} fragment files.".format(len(fragment_files)))
    for file in fragment_files:
        listaBonds = growing_sites(file)
        if firstTime:
            output = input_file(file, user_bond, listaBonds)
            firstTime = False
        else:
            output = input_file(file, user_bond, listaBonds, createFile=False)
    return output.name
