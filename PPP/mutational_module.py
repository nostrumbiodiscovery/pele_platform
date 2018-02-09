import sys

import numpy as np
from prody import writePDB, parsePDB, AtomGroup
from MSM_PELE.PPP.checks_module import CheckMapAndZmatrix
from MSM_PELE.PPP.coordinates_module import GenerateCoordinatesFromZmatrix, ModifyCoords, ModifyCoordinatesPRO
from MSM_PELE.PPP.global_variables import default_mutations_maps
from MSM_PELE.PPP.program_own_classes import ZMATRIX

__author__ = 'jelisa'

# Defining the defaults values for the files thet define the valid maps
# This file should have the format of a python list: ["XXX1", "XXX2", ... , "XXXn"]


def Readmaps(mutation):
    mutation_arbitrary_order = [mutation['ini_resname'], mutation["fin_resname"]]
    mutation_arbitrary_order.sort()
    mutation_maps_key = "-".join(mutation_arbitrary_order)
    try:
        mutation_map = default_mutations_maps[mutation_maps_key]
    except KeyError:
        print "Something is wrong, the mutation you're asking for is not present in the mutations maps."
        sys.exit()
    return mutation_maps_key, mutation_map


def EliminateAtoms(initial_residue, atoms2eliminate):
    """This function gets a residue and a list of atomnames to eliminate from the residue.
    It returns the residue without the atoms."""
    atoms2eliminate = " OR ".join(atoms2eliminate)
    mutated_residue = initial_residue.select('not (name {})'.format(atoms2eliminate)).copy()
    return mutated_residue


def AddAtoms(initial_residue, atoms2add, atomnames_of_2_letters, residue_data, zmatrix, verbose=True):
    """
    This function returns a mutated residue with the atoms added.
    This function gets a residue and a list of atomnames to add to the mutation.
    :rtype : object
    :param initial_residue:
    :param atoms2add: 
    :param atomnames_of_2_letters: 
    :param residue_data: 
    :param zmatrix: 
    :param verbose: 
    """
    if verbose:
        print 'Adding new atoms:'
    new_atoms_elements = []
    for atom in atoms2add:
        if verbose:
            print " # {}".format(atom)
        if atom[:2] in atomnames_of_2_letters:
            new_atoms_elements.append(atom[:2])
        else:
            new_atoms_elements.append(atom[0])

    new_atoms = AtomGroup('New atoms')
    number_of_new_atoms = len(atoms2add)
    new_atoms.setNames(list(atoms2add))
    new_atoms.setElements(new_atoms_elements)
    new_atoms.setResnames([residue_data['fin_resname']] * number_of_new_atoms)
    new_atoms.setResnums([residue_data['resnum']] * number_of_new_atoms)
    new_atoms.setChids([residue_data['chain']] * number_of_new_atoms)
    new_atoms.setAltlocs([''] * number_of_new_atoms)
    new_atoms.setBetas([0] * number_of_new_atoms)
    new_atoms.setIcodes([''] * number_of_new_atoms)
    new_atoms.setOccupancies([1] * number_of_new_atoms)
    new_atoms.setSegnames([''] * number_of_new_atoms)
    new_atoms.setSerials([0] * number_of_new_atoms)
    temporary_coords = np.ones([len(atoms2add), 3])
    new_atoms.setCoords(temporary_coords)
    incomplete_residue = initial_residue + new_atoms
    for atom in new_atoms.iterAtoms():
        atom_new_coordinates = GenerateCoordinatesFromZmatrix(incomplete_residue, [atom.getName()],
                                                                                 zmatrix)
        atom.setCoords(atom_new_coordinates)
        incomplete_residue = initial_residue + new_atoms
    if verbose:
        print 'Done'
    return initial_residue + new_atoms


def ModifyExistingAtoms(initial_residue, atoms2change, atomnames_of_2_letters, initial_atoms_index, final_atoms_index,
                        zmatrix):
    residue2modify = initial_residue.copy()
    print "    * Modifying the atoms:"
    for at_pair in atoms2change:
        ini_atomname = at_pair[initial_atoms_index]
        fin_atomname = at_pair[final_atoms_index]
        print "      # {} -> {}".format(at_pair[initial_atoms_index], at_pair[final_atoms_index])
        element = fin_atomname[:2]
        if element not in atomnames_of_2_letters:
            element = at_pair[final_atoms_index][0]
        atom2modify = residue2modify.select("name {0}".format(ini_atomname))
        atom2modify.setNames([fin_atomname])
        atom2modify.setElements([element])
        if fin_atomname == 'H':
            atom_coords = ModifyCoordinatesPRO(residue2modify, zmatrix, fin_atomname, fin_atomname)
        else:
            atom_coords = GenerateCoordinatesFromZmatrix(residue2modify, [fin_atomname], zmatrix)
        atom2modify.setCoords(atom_coords)
    print "      Done"
    return residue2modify


def ModifyFromToGLY(initial_residue, zmatrix, resname, atoms2change, initial_atoms_index, final_atoms_index, atomnames_of_2_letters):
    modified_residue = initial_residue.copy()
    initial_atoms, final_atoms = [], []
    print "Modifying the atoms:"
    for at_pair in atoms2change:
        initial_atoms.append(at_pair[initial_atoms_index])
        final_atoms.append(at_pair[final_atoms_index])
    for at_pair in atoms2change:
        initial_name = at_pair[initial_atoms_index]
        final_name = at_pair[final_atoms_index]
        print " # {} -> {}".format(initial_name, final_name)
        atom2modify = modified_residue.select("name {}".format(initial_name))
        atom2modify.setNames([final_name])
        if final_name == "H":
            coords = ModifyCoordinatesPRO(modified_residue, zmatrix, final_name, final_name)
        else:
            coords = ModifyCoords(modified_residue, zmatrix, final_name)
        atom2modify.setCoords(coords)
        if final_name in atomnames_of_2_letters:
            pass
        else:
            atom2modify.setElements([final_name[0]])
        atom2modify.setBetas([0])
        atom2modify.setAltlocs([""])
    modified_residue.setResnames([resname] * modified_residue.numAtoms())
    return modified_residue


def AddHfromPRO(initial_residue, zmatrix, mutation):
    new_atom = AtomGroup('Hydrogen')
    new_atom.setNames(['H'])
    new_atom.setElements(['H'])
    new_atom.setResnames([zmatrix.Name])
    new_atom.setResnums([mutation['resnum']])
    new_atom.setChids([mutation['chain']])
    new_atom.setAltlocs([''])
    new_atom.setBetas([0])
    new_atom.setIcodes([''])
    new_atom.setOccupancies([1])
    new_atom.setSegnames([''])
    new_atom.setSerials([0])
    new_atom.setCoords(ModifyCoordinatesPRO(initial_residue, zmatrix, 'H', 'CD'))
    return initial_residue + new_atom

def Mutate(wt_structure, mutation):
    """
    This is the main method that mutates the protein it calls all the other functions needed to do the mutation.
    It reads the mutation map, applies it to the residue, then checks the geometry and coordinates, and finally returns
    the mutated protein.
    """

    # Separate the residue from the rest of the protein, and the rest of the protein is divided in two the part before
    # the residue to mutate, and the part after it. To ease the union of the parts after the mutation process.
    if mutation['chain']:
        initial_residue = wt_structure.select("resname {0[ini_resname]} and chain {0[chain]} and resnum {0[resnum]}".format(mutation)).copy()
        if wt_structure.select("not chain {0[chain]}".format(mutation)) is not None:
            rest_of_the_structure = wt_structure.select("not chain {0[chain]}".format(mutation)).copy()
        else:
            rest_of_the_structure = None
        initial_part_of_the_protein = wt_structure.select("resnum < {0[resnum]} and chain {0[chain]}".format(mutation)).copy()
        final_part_of_the_protein = wt_structure.select("resnum > {0[resnum]} and chain {0[chain]}".format(mutation)).copy()
    else:
        initial_residue = wt_structure.select("resname {0[ini_resname]} and resnum {0[resnum]}".format(mutation)).copy()
        rest_of_the_structure = ''
        initial_part_of_the_protein = wt_structure.select("resnum < {0[resnum]}".format(mutation)).copy()
        final_part_of_the_protein = wt_structure.select("resnum > {0[resnum]}".format(mutation)).copy()


    # This block reads the map and extracts the necessary information
    mutation_maps_key, mutation_map = Readmaps(mutation)
    changing_atoms = mutation_map[1]
    volatile_atoms = mutation_map[2]
    volatile_atoms_behaviour, initial_atoms_index, final_atoms_index = mutation_map[3][
        "-".join([mutation['ini_resname'], mutation['fin_resname']])]

    zmatrix = ZMATRIX(mutation['fin_resname'])
    if zmatrix.Name is None:
        print "The template for the final residue doesn't exist. Check your environment parameters."
        sys.exit()
    CheckMapAndZmatrix(zmatrix.AtomNames, mutation_map, mutation, initial_residue)
    atoms2check_for_clashes = [atom_pair[final_atoms_index] for atom_pair in changing_atoms]

    # This handles the change in the element and atom name.
    atomnames_of_2_letters = ["FE"]  # Modify this list to add more elements that should use 2 letters

    if "GLY" in mutation.values():
        if volatile_atoms_behaviour == "disappear":
            residue2mutate = EliminateAtoms(initial_residue, volatile_atoms)
            mutated_residue = ModifyFromToGLY(residue2mutate, zmatrix, mutation['ini_resname'], changing_atoms,
                                              initial_atoms_index, final_atoms_index, atomnames_of_2_letters)
        elif volatile_atoms_behaviour == "appear":
            residue2mutate = ModifyFromToGLY(initial_residue, zmatrix, mutation['fin_resname'], changing_atoms,
                                             initial_atoms_index, final_atoms_index, atomnames_of_2_letters)
            mutated_residue = AddAtoms(residue2mutate, volatile_atoms, atomnames_of_2_letters, mutation, zmatrix)

        elif volatile_atoms_behaviour == '' and volatile_atoms == []:
            pass
        else:
            print "Something is wrong with the appearing/disappearing value in the file with the maps. Check it!"
            print "The problem is in the entry for the mutation {} it has the value {}. It should be appear or disappear.".format(
                mutation_maps_key, volatile_atoms_behaviour)
            sys.exit()
        if mutation["fin_resname"] == "PRO":
            print 'here'
            atoms2replace = ["CD", "HD2", "HD3"]
            for atom_name in atoms2replace:
                atom = mutated_residue.select('name {}'.format(atom_name))
                coords = GenerateCoordinatesFromZmatrix(mutated_residue, [atom_name], zmatrix)
                atom.setCoords(coords)
        mutated_residue.setResnames([zmatrix.Name] * mutated_residue.numAtoms())
    elif mutation['ini_resname'] == "PRO" and 'H' in volatile_atoms:
        residue2mutate = AddHfromPRO(initial_residue, zmatrix, mutation)
        residue2mutate = ModifyExistingAtoms(residue2mutate, changing_atoms, atomnames_of_2_letters,
                                             initial_atoms_index, final_atoms_index, zmatrix)
        mutated_residue = AddAtoms(residue2mutate, volatile_atoms[1:], atomnames_of_2_letters, mutation, zmatrix)
        mutated_residue.setResnames([zmatrix.Name] * mutated_residue.numAtoms())
    else:
        residue2mutate = ModifyExistingAtoms(initial_residue, changing_atoms, atomnames_of_2_letters,
                                             initial_atoms_index, final_atoms_index, zmatrix)
        residue2mutate.setResnames([mutation['fin_resname']] * initial_residue.numAtoms())
        # Eliminates or Creates the necessary atoms
        if volatile_atoms_behaviour == "disappear":
            mutated_residue = EliminateAtoms(residue2mutate, volatile_atoms)
        elif volatile_atoms_behaviour == "appear":
            mutated_residue = AddAtoms(residue2mutate, volatile_atoms, atomnames_of_2_letters, mutation, zmatrix)
            atoms2check_for_clashes += volatile_atoms
        elif volatile_atoms_behaviour == '' and volatile_atoms == []:
            mutated_residue = residue2mutate
        else:
            print "Something is wrong with the appearing/disappearing value in the file with the maps. Check it!"
            print "The problem is in the entry for the mutation {} it has the value {}. It should be appear or disappear.".format(
                mutation_maps_key, volatile_atoms_behaviour)
            sys.exit()
    if rest_of_the_structure:
        mutated_structure = initial_part_of_the_protein + mutated_residue + final_part_of_the_protein + rest_of_the_structure
    else:
        mutated_structure = initial_part_of_the_protein + mutated_residue + final_part_of_the_protein

    return mutated_structure, zmatrix


def main():
    """
    use like this:
    $> python mutational_module.py pdb_file.pdb 'residue XXX N to YYY'
    """
    initial_structure_path = sys.argv[1]
    mut = sys.argv[2].split()
    mutation = {'ini_resname': mut[1], 'resnum': mut[2], 'fin_resname': mut[4], }
    ff = "force_field parameters"

    initial_structure = parsePDB(initial_structure_path)

    mutated_protein, clashes = Mutate(initial_structure, mutation)
    writePDB('processed.pdb', mutated_protein)
    return mutated_protein, clashes


if __name__ == "__main__":
    mutated_residue, clashes = main()
