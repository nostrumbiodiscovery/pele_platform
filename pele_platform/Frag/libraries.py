import glob
import os


def growing_sites(fragment, user_bond):
    """
    Retrieves all possible growing sites (hydrogens) on the fragment. Takes PDB fragment file as input.
    Output - list of strings represeting sites, e.g. "benzene.pdb C6-H6 C1-H2"
    """
    from rdkit import Chem
    from rdkit.Chem import AllChem
    bonds = []
    mol = Chem.MolFromPDBFile(fragment, removeHs=False)
    print(fragment)
    if mol:
        print("mol", mol)
        print("atoms")
        for a in mol.GetAtoms():
            print(a)
        heavy_atoms = [a for a in mol.GetAtoms() if a.GetSymbol() != "H"]
        for a in heavy_atoms:
            hydrogens = [n for n in a.GetNeighbors() if n.GetSymbol() == "H"]
            at_name = a.GetMonomerInfo().GetName().strip()
            for h in hydrogens:
                h_name = h.GetMonomerInfo().GetName().strip()
                bonds.append("{} {} {}-{}".format(fragment, user_bond, at_name, h_name))
    else:
        print("Couldn't read in molecule from", fragment)
    return bonds


def extract_from_sdf(file_list, path):

    from rdkit import Chem
    output = []

    # setting monomer info to be used for all fragments
    res_info = Chem.AtomPDBResidueInfo()
    res_info.SetResidueName(' GRW')
    res_info.SetChainId('L')
    res_info.SetResidueNumber(1)
    res_info.SetIsHeteroAtom(True)
    elem_count = {}

    for f in file_list:
        mols = Chem.SDMolSupplier(f, removeHs=False)
        for m in mols:
            m.SetProp("_Name", m.GetProp('_Name') + "_converted")
            
            for a in m.GetAtoms():
                n = elem_count.get(a.GetSymbol(), 0)
                n += 1
                elem_count[a.GetSymbol()] = n
                atom_name = " "+a.GetSymbol().strip()+str(n)
                res_info.SetName(atom_name)
                a.SetMonomerInfo(res_info)
            
            # save to PDB file
            file_name = os.path.join(path, "{}.pdb".format(m.GetProp('_Name')))
            writer = Chem.PDBWriter(file_name)
            writer.write(m)
            output.append(file_name)
    return output


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
