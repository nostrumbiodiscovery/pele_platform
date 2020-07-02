import random
import string
import pele_platform.Frag.fragments as fr
import pele_platform.Frag.atoms as at
import pele_platform.Frag.checker as ch
import pele_platform.Errors.custom_errors as ce
import faulthandler; faulthandler.enable()
import time

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
        print("ATOM OF FRAGMENT ATTACHED TO CORE:", atom_fragment)
        print("ATOM OF CORE ATTACHED TO FRAGMENT:", atom_core_idx)
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
        
        for atom in reversed(new_mol.GetMol().GetAtoms()):
            neighbours = atom.GetNeighbors()
            if len(neighbours) == 0: new_mol.RemoveAtom(atom.GetIdx())
        print("FRAGMENT ATOMS", [atom.GetIdx() for atom in new_mol.GetMol().GetAtoms()])
    print("FRAGMENT ATOMS", [atom.GetIdx() for atom in new_mol.GetMol().GetAtoms()])

    # Add missing hydrogen to full ligand and create pdb differently depending on the previous step
    fragment = new_mol.GetMol()
    old_atoms = [atom.GetIdx() for atom in fragment.GetAtoms() if atom.GetAtomicNum() != 1]
    new_atoms = [atom.GetIdx() for atom in ligand.GetAtoms() if
                 atom.GetIdx() not in atoms_core and atom.GetAtomicNum() != 1]
    print("atomsmp", len(old_atoms), len(new_atoms))
    assert len(old_atoms) == len(new_atoms)
    mapping = {new_atom: old_atom for new_atom, old_atom in zip(new_atoms, old_atoms)}
    print("mapping", mapping)
    
    hydrogens_attached_to_fragment_atom = len([atom for atom in ligand.GetAtomWithIdx(atom_fragment).GetNeighbors() if atom.GetAtomicNum() == 1]) + 1
    print("hydrogens_attached_to_fragment_atom", hydrogens_attached_to_fragment_atom)
   
    atom_fragment_mapped = mapping[atom_fragment]
    atom_fragment_valence = fragment.GetAtomWithIdx(atom_fragment_mapped).GetImplicitValence()
    print("atom fragment valence:", atom_fragment_valence, "atomic number:", fragment.GetAtomWithIdx(atom_fragment_mapped).GetAtomicNum(), "mapped atom_fragment id:", atom_fragment_mapped)
    
    # changed July 1st - only this line 
    fragment = Chem.AddHs(fragment, False, True)
    frag_neighbours = fragment.GetAtomWithIdx(atom_fragment_mapped).GetNeighbors()
    lig_neighbours = ligand.GetAtomWithIdx(atom_fragment).GetNeighbors()
    
    print("if condition ", len(frag_neighbours), len(lig_neighbours)) 
    if len(frag_neighbours) < len(lig_neighbours):
        #to_add = hydrogens_attached_to_fragment_atom - atom_fragment_valence
        to_add = len(lig_neighbours) - len(frag_neighbours)
        print("Trying to add H to the fragment...", to_add)
        fragment.GetAtomWithIdx(atom_fragment_mapped).SetNumExplicitHs(to_add)
        t = time.time()
        Chem.MolToPDBFile(fragment, "added_hydrogens_{}.pdb".format(t))
    
    fragment = Chem.AddHs(fragment, False, True)
    t = time.time()
    Chem.MolToPDBFile(fragment, "added_hydrogens2_{}.pdb".format(t)) 
    # checking chirality
    lig_chir = Chem.FindMolChiralCenters(ligand, force=True, includeUnassigned=True)
    print("ligand chirality", ligand, lig_chir)
    frag_chir = Chem.FindMolChiralCenters(fragment, force=True, includeUnassigned=True)
    print("fragment chirality", fragment, frag_chir)

    # check fragment size
    ligand_len = ligand.GetNumHeavyAtoms()
    print("ligand_len", ligand_len)
    frag_len = fragment.GetNumHeavyAtoms()
    print("frag_len", frag_len)
    core_len = ligand_core.GetNumHeavyAtoms()
    print("core_len", core_len)
    if ligand_len - core_len == frag_len:
        if lig_chir and frag_chir:
            for lg, fg in zip(lig_chir, frag_chir):
                #If chirality different or atom number different --> ? as incorrect
                print(mapping)
                print(fg[0])
                print(lg[0])
                if (lg[1] != fg[1] or mapping[lg[0]] != fg[0]):
                    if lg[1] != "?" and fg[1] != "?":
                        print(f"Chirality error {lg[1]} {fg[1]} {mapping[lg[0]]} {fg[0]}")
                        correct = False
                        break
                    else:
                        correct = True
                else:
                    correct = True
        else:
            correct = True
    else:
        print("Atom error")
        correct = False
    return fragment, old_atoms, hydrogen_core, atom_core, atom_fragment, mapping, correct


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
        print("Hydrogen detection failed won't have into account stereochemistry")
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
