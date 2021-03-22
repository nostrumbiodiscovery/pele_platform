import pele_platform.Errors.custom_errors as ce


def check_limit_number_atoms(sdf: str, limit_atoms: int):
    from rdkit import Chem
    """ Function to check input file does not excit certain number of atoms """

    mols = Chem.SDMolSupplier(sdf)
    for mol in mols:
        atoms = mol.GetAtoms()
        natoms = len([atom for atom in atoms])
        if natoms > limit_atoms:
            name_ligand = mol.GetProp("_Name")
            raise ce.LigandSizeExceed(f"Ligand {name_ligand} in {sdf} exceeds size {limit_atoms}")

  
def check_substructure_match(ligand, ligand_core, core_atoms):
    import rdkit.Chem.rdMolAlign as rd
    from rdkit import Chem

    # Initialize ring information
    Chem.GetSSSR(ligand); Chem.GetSSSR(ligand_core)  

    # Align and obtain vocabulary between atoms in core and ligand
    align_obj = rd.GetO3A(ligand, ligand_core)
    atom_map = align_obj.Matches()
    vocabulary = { i2:i1 for i1, i2 in atom_map } 

    # Check the substructure search and change indexes if dont agree with alignment
    for idx_atom_ligand, idx_atom_core in enumerate(core_atoms):
        if idx_atom_ligand not in vocabulary: # If it's a removed hydrogen - pass
            continue
        if vocabulary[idx_atom_ligand] == idx_atom_core: # If indexes are correct - pass
            continue
        else:  # If indexes are incorrect - swap indexes
            core_atoms = list(core_atoms)
            core_atoms[idx_atom_ligand] = vocabulary[idx_atom_ligand]
    return core_atoms
