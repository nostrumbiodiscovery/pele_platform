import glob
import os


def growing_sites(fragment, user_bond):
    """
    Retrieves all possible growing sites (hydrogens) on the fragment. Takes PDB fragment file as input.
    Output - list of strings represeting sites, e.g. "benzene.pdb C6-H6 C1-H2"
    """
    from rdkit import Chem
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


def extract_from_sdf(file_list):

    from rdkit import Chem
    output = []

    for f in file_list:
        mols = Chem.SDMolSupplier(f)
        for m in mols:
            file_name = "{}.pdb".format(m.GetProp('_Name'))
            writer = Chem.PDBWriter(file_name)
            writer.write(m)
            output.extend(file_name)
    return output


def main(user_bond, frag_library):

    # get absolute path to the fragments library
    directory = os.getcwd()
    path = frag_library if os.path.exists(frag_library) else os.path.join(directory, "Libraries", frag_library.strip())
    path = os.path.abspath(path)
    if not os.path.exists(path):
        raise OSError("Cannot find the library.")

    # get fragment files
    fragment_files = []
    extensions = ['*.pdb', '*.sdf', '*.PDB', '*.SDF']
    for e in extensions:
        fragment_files.extend(glob.glob(os.path.join(path, e)))
    print("Retrieved {} fragment files.".format(len(fragment_files)))

    # convert SDF to PDB, if necessary
    sdf_files = [elem for elem in fragment_files if "sdf" in elem.lower()]
    all_files = [elem for elem in fragment_files if "pdb" in elem.lower()]

    if sdf_files:
        all_files.extend(extract_from_sdf(sdf_files))

    # scan fragments for attachment points and write to input.conf
    bond_list = []
    for file in all_files:
        bond_list.extend(growing_sites(file, user_bond))

    output_name = os.path.join(directory, "input.conf")
    with open(output_name, "w+") as conf_file:
        for line in bond_list:
            conf_file.write(line+"\n")

    return output_name
