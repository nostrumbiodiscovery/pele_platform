import sys

from prody import Contacts, calcDistance

from global_processes import FindInitialAndFinalResidues
from global_variables import supported_aminoacids
from program_own_classes import ZMATRIX

__author__ = 'jelisa'

fixed_atoms = ["N", "H", "CA", "HA", "C", "O", "CB"]


def CheckCysteines(structure):
    """
    This function gets one structure and finds the cysteines that form disulfide-bridges and
     the ones that should be charged and returns their residue numbers and chain.

    :param structure: Prody Atomgroup object
    :return: two different lists
    """
    crosslinked_cys = []
    charged_cys = []
    if structure.select('resname CYS') is not None:
        for cysteine in structure.select('resname CYS').copy().iterResidues():
            if "HG" in cysteine.getNames():
                continue
            key = "{}_{}".format(cysteine.getResnum(), cysteine.getChid())
            if key in crosslinked_cys:
                continue
            elif key in charged_cys:
                print "CYS/CYT classification WTF.0??\n cystein {} {}".format(cysteine.getResnum(),
                                                                              cysteine.getChid())
                continue
            sg = cysteine.select('name SG')
            if sg is None:
                continue
            selection_string = "(within 2.1 of atom) and (not (resnum {} and name SG))".format(cysteine.getResnum())
            neighbours = structure.select(selection_string, atom=sg)
            if neighbours is None:
                print "CYS/CYT classification WTF.1??"
            elif "SG" in neighbours.getNames():
                crosslinked_cys.append(key)
                key2 = "{}_{}".format(neighbours.select('name SG').getResnums()[0],
                                      neighbours.select('name SG').getChids()[0])
                crosslinked_cys.append(key2)
            else:
                charged_cys.append(key)
    return crosslinked_cys, charged_cys


def CheckStructure(initial_structure, gaps={}, no_gaps={}, remove_missing_ter=False, debug=False):
    residues2fix = {}
    crosslinked_cysteines, charged_cysteines = CheckCysteines(initial_structure)
    residues2remove = {}
    missing_residues = []
    for chain in initial_structure.iterChains():
        if chain.getChid() in gaps.keys():
            gaps_residues = gaps[chain.getChid()]
        else:
            gaps_residues = []
        if chain.getChid() in no_gaps.keys():
            no_gaps_residues = no_gaps[chain.getChid()]
        else:
            no_gaps_residues = []
        initial_residue, final_residue = FindInitialAndFinalResidues(chain)
        residues2ignore = []
        for residue in chain.iterResidues():
            resname = residue.getResname().strip()
            resnum = residue.getResnum()
	    reschain = residue.getChid()
            residue_atomnames = list(residue.getNames())
            if resname[:2] == "CU":
                resname = "CU"
            if resnum == initial_residue and resname in supported_aminoacids:
                zmatrix = ZMATRIX(resname + 'B')
            elif resnum == final_residue:
                zmatrix = ZMATRIX(resname + "E")
            elif resname == "NMA":
                zmatrix = ZMATRIX(resname + "E")
            elif resname == "ACE":
                zmatrix = ZMATRIX(resname + "B")
            else:
                try:
                    zmatrix = ZMATRIX(resname)
                except:
                    print "UP {}".format(resname)
            if zmatrix.Name is None:
                print "  The residue {} {} doesn't have a template, so it won't be checked.".format(resname, resnum)
                missing_residues.append([resname, reschain])
                continue
            residue_atomnames.sort()
            if sorted(zmatrix.AtomNames) != residue_atomnames:
                missing_atoms = [atom_name for atom_name in zmatrix.AtomNames if atom_name not in residue_atomnames]
                over_atoms = [atom_name for atom_name in residue_atomnames if atom_name not in zmatrix.AtomNames]
                if over_atoms:
                    if residue.getResnum() in gaps_residues or residue.getResnum() in [initial_residue, final_residue]:
                        if "HXT" in over_atoms:
                            over_atoms.pop(over_atoms.index('HXT'))
                        else:
                            if "H1" in over_atoms:
                                over_atoms.pop(over_atoms.index('H1'))
                            if "H2" in over_atoms:
                                over_atoms.pop(over_atoms.index('H2'))
                    else:
                        print "   The residue {} {} has more atoms than the zmatrix.".format(resname, resnum)
                    if debug:
                        print 'over_atoms: ', over_atoms
                        print "ZMATRIX: {}\n zmat atoms:{}".format(zmatrix.Name, sorted(zmatrix.AtomNames))
                        print " resi atoms:", sorted(residue_atomnames)
                if missing_atoms:
                    if "CA" in missing_atoms or "N" in missing_atoms or "O" in missing_atoms or "C" in missing_atoms:
                        print "  The residue {} {} is missing one or more heavy atoms in the backbone. " \
                              "It cannot be fixed.".format(residue.getResname(), residue.getResnum(), residue.getChid())
                        if remove_missing_ter:
                            if residue.getResnum() in [initial_residue, final_residue]:
                                print " INFO: In structure {} the residue {} {} in chain {} will be eliminated".format(
                                    initial_structure.getTitle(), residue.getResname(),
                                    residue.getResnum(), residue.getChid())
                                try:
                                    residues2remove[residue.getChid()]
                                except KeyError:
                                    residues2remove[residue.getChid()] = [residue.getResnum()]
                                else:
                                    residues2remove[residue.getChid()].append(residue.getResnum())
                                if resnum == initial_residue:
                                    next_res_resnum = resnum + 1
                                    next_residue = chain.getResidue(next_res_resnum)
                                    next_res_resname = next_residue.getResname()
                                    if next_res_resname != "PRO":
                                        atom = next_residue.getAtom("H")
                                        if atom is None:
                                            atoms2add = {"H1", "H2", "H3"}
                                        else:
                                            residues2ignore.append(next_res_resnum)
                                            atom.setName("H1")
                                            atoms2add = {"H2", "H3"}
                                    else:
                                        atom = residue.getAtom("H")
                                        if atom is None:
                                            atoms2add = {"H1", "H2"}
                                        else:
                                            atom.setName("H1")
                                            atoms2add = {"H2"}
                                    residues2fix["{} {} {}".format(next_res_resname, next_res_resnum,
                                                                   residue.getChid())] = atoms2add
                                elif resnum == final_residue:
                                    prev_res_resnum = resnum - 1
                                    prev_residue = chain.getResidue(prev_res_resnum)
                                    if prev_residue is None:
                                        print "There's no residue with number {} in chain {}".format(prev_res_resnum,
                                                                                                     chain.getChid())
                                    prev_res_resname = prev_residue.getResnames()[0]
                                    atoms2add = {"OXT"}
                                    residues2fix["{} {} {}".format(prev_res_resname, prev_res_resnum,
                                                                   residue.getChid())] = atoms2add
                        else:
                            print "  The program will be interrupted.".format(resname, resnum)
                            print "     This are the missing atoms: {}".format(missing_atoms)
                            sys.exit()
                    else:
                        if "H" in missing_atoms and resnum in gaps_residues and resnum not in residues2ignore:
                            maestro_terminal_H = ["H1", "H2", "H11", "H22"]
                            atomname_to_use = [atom_name for atom_name in maestro_terminal_H
                                               if atom_name in residue.getNames()][0]
                            if atomname_to_use:
                                atom = residue.getAtom(atomname_to_use)
                                atom.setName('H')
                                missing_atoms.pop(missing_atoms.index('H'))

                        key = "{}_{}".format(resnum, residue.getChid())
                        if resname == "CYS":
                            if key in crosslinked_cysteines + charged_cysteines:
                                if key in charged_cysteines:
                                    residue.setResname('CYT')
                                if missing_atoms == ['HG']:
                                    continue
                                elif "HG" in missing_atoms:
                                    missing_atoms.pop(missing_atoms.index('HG'))
                        for atom_name in missing_atoms:
                            if atom_name[0] != "H" and atom_name[0].isalpha():
                                print "  The residue {} {} is missing the heavy atom {}.\n" \
                                      "  All the atoms depending o this atom will be placed according to the zmatrix.".format(
                                        resname, resnum, atom_name)
                                atoms2add = zmatrix.GetAllChildrenAtoms(atom_name)
                                residues2fix["{} {} {}".format(resname, resnum, residue.getChid())] = set(atoms2add).union(
                                    set(missing_atoms))
                            else:
                                if resnum not in residues2ignore:
                                    residues2fix["{} {} {}".format(resname, resnum, residue.getChid())] = set(
                                        missing_atoms)
    return residues2fix, residues2remove, missing_residues


def CheckMapAndZmatrix(zmap_atoms, mutation_map, mutation, residue):
    """This function checks that the number of atoms in the residue to mutate agree with the number of atoms in the map.
    It also checks the agreement between the mutation map and the zmatrix"""

    common_atoms = mutation_map[0]
    changing_atoms = mutation_map[1]
    volatile_atoms = mutation_map[2]
    volatile_atoms_behaviour, initial_atoms_index, final_atoms_index = mutation_map[3][
        "-".join([mutation['ini_resname'], mutation['fin_resname']])]

    # A small check to see that the map and the residue agree in the number of atoms.
    map_number_of_atoms = len(common_atoms) + len(changing_atoms) + len(volatile_atoms)
    residue2mutate_number_of_atoms = residue.numAtoms()
    map_atoms = common_atoms + [atom_pair[final_atoms_index] for atom_pair in changing_atoms]

    initial_residue_map_atoms = common_atoms + [atom_pair[initial_atoms_index] for atom_pair in changing_atoms]
    if volatile_atoms_behaviour == 'disappear':
        initial_residue_map_atoms += volatile_atoms
    else:
        map_atoms += volatile_atoms
    if map_number_of_atoms != residue2mutate_number_of_atoms and map_number_of_atoms != residue2mutate_number_of_atoms + len(
            volatile_atoms):
        print "Check the number of atoms of the structure because it doesn't match with the map."
        print "Remember that mutations at the beginning or ending of the protein are not supported yet."
        print "residue atoms: {}".format(set(residue.getNames()))
        print "initial residue map atoms: {}".format(initial_residue_map_atoms)
        sys.exit()

    if volatile_atoms_behaviour == 'appear':
        map_atoms += volatile_atoms
    if set(map_atoms) != set(zmap_atoms):
        print "The zmatrix for the desired aminoacid and the map for the mutation don't agree. Check Them."
        print "map atoms: {}".format(set(map_atoms))
        print "zmatrix atoms: {}".format(set(zmap_atoms))
        sys.exit()


def CheckMutation(wt_structure, mutation):
    """
     This function checks whether the specified mutation is possible to make or not. And specifies the reason why it
    doesn't work.
     It may not be possible to make it due to the fact that the specified residue isn't present in the protein, or
    the final aminoacid to generate maybe isn't supported by the program.
    This function has as input the:
        - wt_structure: the structure to mutate
        - mutation: the mutation to do with the format ["original_resname", "resnum", "final_resname"]
        - valid_resnames_file: the path to the file containing a list with the supported mutations
    """

    if mutation['fin_resname'] not in supported_aminoacids:
        print "The desired mutation isn't possible because the desired aa. is not supported."
        print supported_aminoacids
        return False

    testing_the_resnum = "resnum {}".format(mutation['resnum'])
    testing_the_resname = "resname {}".format(mutation['ini_resname'])

    if wt_structure.select(testing_the_resnum) is None:
        print "The protein doesn't have a residue numbered {}".format(mutation['resnum'])
        return False
    elif wt_structure.select(testing_the_resname) is None:
        print "The protein doesn't have a residue named {}".format(mutation['ini_resname'])
        # print wt_structure.getResnames()
        return False
    elif mutation['chain']:
        testing_the_chain = "chain `{}` ".format(mutation['chain'])
        if wt_structure.select(testing_the_chain) is None:
            print "The protein doesn't have a chain named {}".format(mutation['chain'])
        elif wt_structure.select(testing_the_resname + ' and ' + testing_the_resnum + ' and ' + testing_the_chain) is None:
            print "The protein doesn't have the residue {0[ini_resname]} {0[chain]} {0[resnum]}".format(mutation)
            return False
    else:
        return True


def CheckClashes(mutated_protein, mutation, zmatrix, initial_residue,
                 final_residue, overlapfactor=0.7):
    """
    :param initial_residue:
    :param mutation:
    :rtype : dictionary
    :param mutated_protein: prody AtomGroup
    :param zmatrix: ZMATRIX
    :param overlapfactor: float

    A function to check the clashes between a given residue and the rest of the
    structure. It needs the number of the mutated residue the structure, the
    zmatrix of the mutated residue and the overlapfactor.

    The mutated protein is an Atomgroup object from prody that has the structure
    to check for clashes.

    The num_mutated_residue is the number of the mutated residue to check for
    clashes with the rest of the structure.

    The zmatrix contains the information extracted from the template of PELE.

    The overlapfactor by default is 0.7 but it can be set by the user to any
    number between 0 and 1. If it's one the atoms can be completely superimposed
    if the number is 0 they can't superimpose at all.
    """
    max_vdw = max(zmatrix.VdWRadius)
    contacts_object = Contacts(mutated_protein)
    if mutation["chain"]:
        mutated_residue = mutated_protein.select("resnum {0[resnum]} and chain {0[chain]} and not hydrogen".format(mutation))
    else:
        mutated_residue = mutated_protein.select("resnum {0[resnum]} and not hydrogen".format(mutation))
    mutated_residue_indices = mutated_residue.getIndices()
    clashes2check = {}
    for at in mutated_residue.iterAtoms():
        atom_name = at.getName()
        if atom_name in fixed_atoms:
            continue
        try:
            at_vdw = zmatrix.VdWRadius[zmatrix.AtomNames.index(atom_name)]
        except ValueError:
            print "The atom {} in the mutated residue doesn't agree with the Zmatrix.".format(atom_name)
            print mutated_residue.getNames()
            print zmatrix.AtomNames
            raise ValueError
        dist = (at_vdw + max_vdw) * overlapfactor
        at_contacts = contacts_object.select(dist, at.getCoords())
        if at_contacts is None:
            continue
        possible_clashes = [clashing_atom.getIndex() for clashing_atom in at_contacts
                            if clashing_atom.getIndex() not in mutated_residue_indices]
        if possible_clashes:
            clashes2check[atom_name] = possible_clashes

    if mutated_residue.getResnames()[0] == "PRO":
        try:
            cd_atom_clashes = clashes2check["CD"]
        except KeyError:
            pass
        else:
            n_atom_index = mutated_residue.select("name N").getIndices()[0]
            try:
                cd_atom_clashes.index(n_atom_index)
            except ValueError:
                pass
            else:
                cd_atom_clashes.pop(cd_atom_clashes.index(n_atom_index))
    real_clashes = {}
    for key, indices in clashes2check.iteritems():
        real_clashing_indices = []
        for ind in indices:
            atom2check = mutated_protein.select('index {}'.format(ind))
            atom2check_resname = atom2check.getResnames()[0]
            atom2check_resnum = atom2check.getResnums()[0]
            if atom2check_resnum == initial_residue and atom2check_resname in supported_aminoacids:
                atom2check_zmatrix = ZMATRIX(atom2check_resname + 'B')
            elif atom2check_resnum == final_residue:
                atom2check_zmatrix = ZMATRIX(atom2check_resname + 'E')
            else:
                atom2check_zmatrix = ZMATRIX(atom2check_resname)
            if atom2check_zmatrix.Name is None:
                print "Cannot load the zmatrix for the residue {} {}.\n" \
                      "It won't be checked for clashes.".format(atom2check_resname, atom2check_resnum)
                continue
            atom2check_vdw_radius = atom2check_zmatrix.VdWRadius[
                atom2check_zmatrix.AtomNames.index(atom2check.getNames()[0])]
            key_vwd_radius = zmatrix.VdWRadius[zmatrix.AtomNames.index(key)]
            min_distance = (atom2check_vdw_radius + key_vwd_radius) * overlapfactor
            real_distance = calcDistance(mutated_residue.select('name {}'.format(key)), atom2check)
            if real_distance < min_distance:
                real_clashing_indices.append(ind)
        if real_clashing_indices:
            real_clashes[key] = real_clashing_indices

    return real_clashes


def CheckforGaps(structure):
    gaps = {}
    not_gaps = {}
    for chain in structure.iterChains():
        # initial, final = FindInitialAndFinalResidues(chain)
        previous_residue_number = None
        for residue in chain.iterResidues():
            chain_id = chain.getChid()
            # print residue.getResnum() < previous_residue_number
            if previous_residue_number is not None:
                if residue.getResnum() > previous_residue_number + 1:
                    current_residue_N = residue.getAtom("N")
                    previous_residue = chain.getResidue(previous_residue_number)
                    previous_residue_C = previous_residue.getAtom('C')
                    if residue.getResname() in supported_aminoacids:
                        if current_residue_N is not None and previous_residue_C is not None:
                            distance = calcDistance(current_residue_N, previous_residue_C)
                            if distance < 1.5:
                                try:
                                    not_gaps[chain_id]
                                except KeyError:
                                    not_gaps[chain_id] = []
                                not_gaps[chain_id].extend([previous_residue_number, residue.getResnum()])
                            else:
                                try:
                                    gaps[chain_id]
                                except KeyError:
                                    gaps[chain_id] = []
                                gaps[chain_id].extend([previous_residue_number, residue.getResnum()])
                        elif current_residue_N is None:
                            print "There's a problem with residue {} {} {} it" \
                                  " doesn't have the N atom".format(residue.getResname(),
                                                                    residue.getResnum(),
                                                                    residue.getChid())
                        elif previous_residue_C is None:
                            print "There's a problem with residue {} {} {} it doesn't have the C atom".format(
                                previous_residue.getResname(),
                                previous_residue.getResnum(),
                                previous_residue.getChid())
                elif residue.getResnum() == previous_residue_number:
                    if residue.getIcode() == '':
                        if residue.hetero is not None:
                            residue.setResnum(previous_residue_number + 1)
                        else:
                            print 'WARNING: There are two residues with the same resnum in the same chain.' \
                                  'The residue {} {} {} '.format(residue.getChid(),
                                                                 residue.getResnum(),
                                                                 residue.getResname())
                elif residue.getResnum() < previous_residue_number + 1:
                    if residue.hetero is not None:
                        residue.setResnum(previous_residue_number + 1)
                    else:
                        print residue.getResnum(), previous_residue_number
                        print "WARNING! The next residue has a lower resnum than the previous one in the chain??"
                        print "CHAIN: {} RESIDUE: {}  RESNUM {}".format(residue.getChid(),
                                                                        residue.getResnum(),
                                                                        residue.getResname())
                        print "previous: RESIDUE {} RESNUM  {}".format(previous_residue.getResname(),
                                                                       previous_residue.getResnum())
                        return None, None
                else:
                    pass
            previous_residue_number = residue.getResnum()
    return gaps, not_gaps

