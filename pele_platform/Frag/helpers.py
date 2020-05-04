import os
import pele_platform.Frag.fragments as fr
import pele_platform.Frag.atoms as at
import pele_platform.Frag.checker as ch
import pele_platform.Errors.custom_errors as ce



def _search_core_fragment_linker(ligand, ligand_core, result=0, check_simmetry=False):
    """ Given mol1 and mol2 return the linker atoms"""
    substructure_results = ligand.GetSubstructMatches(ligand_core)
    try:
        core_atoms = substructure_results[result]
    except IndexError:
        raise IndexError("Make sure core from pdb and full fragment are the same. Be carefull \
that either core and fragment have corectly define aromatic bonds. Also, all fragments must have different molecule name")
    # Sometime substructure search mess up with symettry. Check that!
    if check_simmetry:
        core_atoms = ch.chec_substructure_match(ligand, ligand_core, core_atoms) 
    for atom in ligand.GetAtoms():
        if atom.GetIdx() in core_atoms or atom.GetAtomicNum() == 1:
            continue
        else:
            atoms_bonded = atom.GetNeighbors()
            for neighbour in atoms_bonded:
                if neighbour.GetIdx()  in core_atoms:
                    return core_atoms.index(neighbour.GetIdx()), core_atoms, substructure_results


def _build_fragment_from_complex(complex, residue, ligand, ligand_core, result=0, substructure=True, simmetry=False):
    from rdkit import Chem
    import rdkit.Chem.rdmolops as rd
    import rdkit.Chem.rdchem as rc
    import rdkit.Chem.AllChem as rp

    #Retrieve atom core linking fragment
    try:
        atom_core_idx, atoms_core, _ = _search_core_fragment_linker(ligand, ligand_core, result, simmetry)
    except TypeError:
        raise ce.SameMolecule("Core and ligand are the exact same molecule. Check your inputs")
    atom_core = at.Atom(ligand_core, atom_core_idx)
    mol = Chem.MolFromPDBFile(complex, removeHs=False)

    #Retrieve hydrogen core attach to linking atom
    original = rd.SplitMolByPDBResidues(mol)[residue]
    hydrogen_core_idx = [atom.GetIdx() for atom in original.GetAtomWithIdx(atom_core_idx).GetNeighbors() if atom.GetAtomicNum() == 1][0]
    hydrogen_core = at.Atom(original, hydrogen_core_idx)

    # Delete core for full ligand with substructure 
    # and if it fails manually
    #Chem.MolToPDBFile(ligand, "int0.pdb")
    if substructure:
        fragment = rd.DeleteSubstructs(ligand, ligand_core)
        new_mol = rc.EditableMol(fragment)
        for atom in reversed(fragment.GetAtoms()):
            neighbours = atom.GetNeighbors()
            if len(neighbours) == 0: new_mol.RemoveAtom(atom.GetIdx())
    else:
        new_mol = rc.EditableMol(ligand)
        for atom in atoms_core:
            new_mol.RemoveAtom(atom)
        #Chem.MolToPDBFile(new_mol.GetMol(), "int1.pdb")
        for atom in reversed(new_mol.GetMol().GetAtoms()): 
            neighbours = atom.GetNeighbors()
            if len(neighbours) == 0: new_mol.RemoveAtom(atom.GetIdx())


    #Add missing hydrogen to full ligand and create pdb differently
    #depending on the previous step
    fragment = new_mol.GetMol()
    #Chem.MolToPDBFile(fragment, "int2.pdb")
    old_atoms = [atom.GetIdx() for atom in fragment.GetAtoms()]
    if substructure:
        fragment = Chem.AddHs(fragment, False, True)
    else:
        #res = rp.EmbedMolecule(fragment)
        fragment = Chem.AddHs(fragment, False, True)
    return fragment, old_atoms, hydrogen_core, atom_core

def _retrieve_fragment(fragment, old_atoms, atom_core, hydrogen_core): 
    from rdkit import Chem
    import rdkit.Chem.rdmolops as rd
    import rdkit.Chem.rdchem as rc
    # Get new added hydrogen
    try:
        added_hydrogen_idx = [atom.GetIdx() for atom in fragment.GetAtoms() if atom.GetIdx() not in old_atoms][0]
        no_hydrogens = False
    except IndexError:
        print("Hydrogen detection failed won't have into account steriochemistry")
        added_hydrogen_idx = 0
        no_hydrogens = True
    # Get fragment atom attached to newly added hydrogen
    atom_fragment_idx = fragment.GetAtomWithIdx(added_hydrogen_idx).GetNeighbors()[0].GetIdx()
    # Save fragment
    fragment_filename = fragment.GetProp("_Name")+".pdb"
    Chem.MolToPDBFile(fragment, fragment_filename)
    fragment = Chem.MolFromPDBFile(fragment_filename, removeHs=False)
    added_hydrogen = at.Atom(fragment, added_hydrogen_idx)
    atom_fragment_attach_to_hydrogen = at.Atom(fragment, atom_fragment_idx)
    # Build fragment object
    fragment = fr.Fragment(fragment_filename, atom_fragment_attach_to_hydrogen, 
        added_hydrogen, atom_core, hydrogen_core, no_hydrogens)
    return fragment
