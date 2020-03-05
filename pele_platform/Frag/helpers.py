



def search_core_fragment_linker(ligand, ligand_core, result=0):
    """ Given mol1 and mol2 return the linker atoms"""
    substructure_results = ligand.GetSubstructMatches(ligand_core)
    core_atoms = substructure_results[result]
    for atom in ligand.GetAtoms():
        if atom.GetIdx() in core_atoms:
            continue
        else:
            atoms_bonded = atom.GetNeighbors()
            for neighbour in atoms_bonded:
                if neighbour.GetIdx()  in core_atoms:
                    return core_atoms.index(neighbour.GetIdx()), core_atoms, substructure_results
