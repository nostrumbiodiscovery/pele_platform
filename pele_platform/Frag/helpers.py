import random
import string
import pele_platform.Frag.fragments as fr
import pele_platform.Frag.atoms as at
import pele_platform.Frag.checker as ch
import pele_platform.Errors.custom_errors as ce


def _search_core_fragment_linker(ligand, ligand_core, result=0, check_symmetry=False):
    """ Given mol1 and mol2 return the linker atoms"""
    substructure_results = ligand.GetSubstructMatches(ligand_core)

    try:
        print("Running substucture search")
        core_atoms = substructure_results[result]
        print("RESULT SUBSTRUCTURE", core_atoms)
    except IndexError:
        raise IndexError("Make sure core from pdb and full fragment are the same. Be careful \
that either core and fragment have corectly define aromatic bonds. Also, all fragments must have different molecule name")
    # Sometime substructure search mess up with symmetry. Check that!
    if check_symmetry:
        core_atoms = ch.check_substructure_match(ligand, ligand_core, core_atoms)
    for atom in ligand.GetAtoms():
        if atom.GetIdx() in core_atoms or atom.GetAtomicNum() == 1:
            continue
        else:
            atoms_bonded = atom.GetNeighbors()
            for neighbour in atoms_bonded:
                if neighbour.GetIdx() in core_atoms:
                    return core_atoms.index(neighbour.GetIdx()), core_atoms, atom.GetIdx()


def _build_fragment_from_complex(complex, residue, ligand, ligand_core, result=0, substructure=True, simmetry=False):
    from rdkit import Chem
    import rdkit.Chem.rdmolops as rd
    import rdkit.Chem.rdchem as rc

    # Retrieve atom core linking fragment
    try:
        atom_core_idx, atoms_core, atom_fragment = _search_core_fragment_linker(ligand, ligand_core, result, simmetry)
        print("ATOM OF FRAGMENT ATTACH TO CORE:", atom_fragment)
        print("ATOM OF CORE ATTACH TO FRAG:", atom_core_idx)
    except TypeError:
        raise ce.SameMolecule("Core and ligand are the exact same molecule. Check your inputs")
    atom_core = at.Atom(ligand_core, atom_core_idx)
    mol = Chem.MolFromPDBFile(complex, removeHs=False)

    # Retrieve hydrogen core attach to linking atom
    original = rd.SplitMolByPDBResidues(mol)[residue]
    hydrogen_core_idx = \
    [atom.GetIdx() for atom in original.GetAtomWithIdx(atom_core_idx).GetNeighbors() if atom.GetAtomicNum() == 1][0]
    hydrogen_core = at.Atom(original, hydrogen_core_idx)

    # Delete core for full ligand with substructure and if it fails manually
    print("MOL ATOMS", [atom.GetIdx() for atom in ligand.GetAtoms()])
    if substructure:
        print("Remove fragment from full ligand 1")
        fragment = rd.DeleteSubstructs(ligand, ligand_core)
        atoms = [i.GetIdx() for i in fragment.GetAtoms()]
        print("Fragment atoms", atoms)
        new_mol = rc.EditableMol(fragment)
        for atom in reversed(fragment.GetAtoms()):
            neighbours = atom.GetNeighbors()
            if len(neighbours) == 0: new_mol.RemoveAtom(atom.GetIdx())
    else:
        print("Remove fragment from full ligand 2")
        new_mol = rc.EditableMol(ligand)
        Chem.MolToPDBFile(new_mol.GetMol(), "int1.0.pdb")
        print(atoms_core)
        for atom in reversed(atoms_core):
            print("AAA", atom)
            new_mol.RemoveAtom(atom)
        print("B")
        print("C")
        for atom in reversed(new_mol.GetMol().GetAtoms()):
            neighbours = atom.GetNeighbors()
            if len(neighbours) == 0: new_mol.RemoveAtom(atom.GetIdx())
        print("FRAGMENT ATOMS", [atom.GetIdx() for atom in new_mol.GetMol().GetAtoms()])

    # Add missing hydrogen to full ligand and create pdb differently depending on the previous step
    fragment = new_mol.GetMol()
    old_atoms = [atom.GetIdx() for atom in fragment.GetAtoms()]
    new_atoms = [atom.GetIdx() for atom in ligand.GetAtoms() if
                 atom.GetIdx() not in atoms_core and atom.GetAtomicNum() != 1]
    mapping = {new_atom: old_atom for new_atom, old_atom in zip(new_atoms, old_atoms)}
    try:
        Chem.MolToPDBFile(fragment, "int2.pdb")
    except Exception:
        hydrogens_attached_to_fragment_atom = len(
            [atom for atom in ligand.GetAtomWithIdx(atom_fragment).GetNeighbors() if atom.GetAtomicNum() == 1]) + 1
        fragment.GetAtomWithIdx(mapping[atom_fragment]).SetNumExplicitHs(hydrogens_attached_to_fragment_atom)
    if substructure:
        fragment = Chem.AddHs(fragment, False, True)
    else:
        fragment = Chem.AddHs(fragment, False, True)
    return fragment, old_atoms, hydrogen_core, atom_core, atom_fragment, mapping


def random_string(string_length=8):
    letters = string.ascii_lowercase
    return ''.join(random.choice(letters) for i in range(string_length))


def _retrieve_fragment(fragment, old_atoms, atom_core, hydrogen_core, atom_fragment, mapping):
    from rdkit import Chem

    # Get new added hydrogen
    try:
        added_hydrogen_idx = [atom.GetIdx() for atom in fragment.GetAtoms() if atom.GetIdx() not in old_atoms][0]
        no_hydrogens = False
        # Get fragment atom attached to newly added hydrogen
        atom_fragment_idx = fragment.GetAtomWithIdx(added_hydrogen_idx).GetNeighbors()[0].GetIdx()
    except IndexError:
        print("Hydrogen detection failed won't have into account steriochemistry")
        added_hydrogen_idx = 0
        no_hydrogens = True
        atom_fragment_idx = mapping[atom_fragment]
        print("Final atom fragment", atom_fragment_idx)

    # Save fragment
    string = random_string(8)
    fragment_filename = fragment.GetProp("_Name") + string + ".pdb"
    print("file fragment", fragment_filename)
    Chem.MolToPDBFile(fragment, fragment_filename)
    fragment = Chem.MolFromPDBFile(fragment_filename, removeHs=False)
    added_hydrogen = at.Atom(fragment, added_hydrogen_idx)
    atom_fragment_attach_to_hydrogen = at.Atom(fragment, atom_fragment_idx)

    # Build fragment object
    fragment = fr.Fragment(fragment_filename, atom_fragment_attach_to_hydrogen, added_hydrogen, atom_core,
                           hydrogen_core, no_hydrogens)
    return fragment
