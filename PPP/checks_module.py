import sys

from prody import Contacts, calcDistance, calcAngle

from global_processes import FindInitialAndFinalResidues
from global_variables import supported_aminoacids, supported_metals, coordination_geometries
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
            key = "{}_{}".format(cysteine.getChid(), cysteine.getResnum())
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
                key2 = "{}_{}".format(neighbours.select('name SG').getChids()[0],
                                      neighbours.select('name SG').getResnums()[0])
                crosslinked_cys.append(key2)

            else:
                charged_cys.append(key)
    if crosslinked_cys:
        writing_list = ["     * CYS_{} - CYS_{}".format(x, y) for x,y in zip(crosslinked_cys[0::2], crosslinked_cys[1::2])]
        print "  * The following cysteins will be considered as cross-linked:\n{0}".format("\n".join(writing_list))
    if charged_cys:
        writing_list = ["     * CYS_{}".format(k) for k in charged_cys]
        print "  * The following cysteins will be considered as cross-linked:\n{0}".format("\n".join(writing_list))
    return crosslinked_cys, charged_cys


# def CheckTetrahedricConformation(angles):
#     right_conformation = False
#     distorted_angles = []
#     ok_angles = []
#     accepted_angles = coordination_geometries['tetrahedric']
#     for angle in angles:
#         if angle[3] < accepted_angles[0] - 15 or angle[3] > accepted_angles[0] + 15:
#             distorted_angles.append(angle)
#         else:
#             ok_angles.append(angle)
#     if len(ok_angles) > len(distorted_angles):
#         if distorted_angles:
#             print "       * The following angle(s) are distorted for this geometry:\n{}".format(
#                 "       * {} {} {}".format(angle[0].getName(), angle[0].getResnum(), angle[0].getChid()))
#         right_conformation = True
#     return  right_conformation


def CheckConformation(angles, conformation):
    """
    A function to check that the angles are within the right range for a given conformation
    :param angles: a list containing the three atoms involved in the angle (prody instances) and the
    angle in the position 3
    :param conformation: a string containing the name of the conformation to check, it should be one defined in
    global_variables.coordination_geometries
    :return: a boolean indicating whether the conformation is right or not
    """
    right_conformation = False
    distorted_angles = []
    ok_angles = []
    accepted_angles, needed_atoms = coordination_geometries[conformation]
    if conformation not in coordination_geometries.keys():
        print "      * The desired coordination isn't implemented."
    else:
        if len(angles) < needed_atoms:
            print "       * WARNING: The metal is missing atoms to coordinated with"

    for angle in angles:
        for a_angle in accepted_angles:
            if a_angle - 15 <= angle[3] <= a_angle + 15:
                ok_angles.append(angle)
            else:
                distorted_angles.append(angle)
    if len(ok_angles) > len(distorted_angles):
        if distorted_angles and len(angles) != needed_atoms:
            print "       * The angles match the {} coordination".format(conformation)
            # print [[[x.getName(), x.getResnum(), x.getChid()] for x in angle[]] for angle in distorted_angles]
            print "       * The angle(s) formed by the following atoms are distorted for this geometry:"#\n{}".format(
            for angle in distorted_angles:
                print "{} {} {} {}".format("        * ",
                                           "{}-{}-{}".format(angle[0].getName(), angle[0].getResnum(), angle[0].getChid(),),
                                           "{}-{}-{}".format(angle[1].getNames()[0], angle[1].getResnums()[0], angle[1].getChid()),
                                           "{}-{}-{}".format(angle[2].getName(), angle[2].getResnum(), angle[2].getChid()))
        else:
            print "       * The angles match perfectly the {} coordination".format(conformation)
        right_conformation = True
    return right_conformation


def CheckMetalsCoordination(structure):
    """
    A function to detect the metal atoms that could be coordinated with the protein.
    :param structure:
    :return:
    """
    selection_pattern = "(within 3 of metal) and (not resnum {}) and (not hydrogen) and (not carbon)"
    coordinated_metals = {}
    for metal in supported_metals:
        if structure.select('resname {}'.format(metal)) is not None:
            print "  * Checking the metals that can be coordinated. (" \
                  "a constraint should be used if they're really coordinated)"
            for metal_res in structure.select('resname {}'.format(metal)).copy().iterResidues():
                coordinated_atoms_list = []
                if metal_res.numAtoms() != 1:
                    continue
                coordinated_atoms = structure.select(selection_pattern.format(metal_res.getResnum()),
                                                     metal=metal_res)
                if coordinated_atoms is None:
                    print "     * The metal atom {} isn't coordinated with the protein. Are you sure it's necessary?"
                else:
                    prev_atom = None
                    for at in coordinated_atoms.iterAtoms():
                        if prev_atom is None:
                            coordinated_atoms_list.append(at)
                            prev_atom = at
                        else:
                            if at.getResnum() == prev_atom.getResnum() and \
                                    (at.getIndex() == prev_atom.getIndex() + 1 or at.getResname() in ['ASP', 'GLU']):
                                if calcDistance(at, metal_res) > calcDistance(prev_atom, metal_res):
                                    prev_atom = None
                                    continue
                                else:
                                    coordinated_atoms_list.pop(-1)
                            coordinated_atoms_list.append(at)
                            prev_atom = at
                    coordinated_metals[metal_res] = coordinated_atoms_list#, 'angles': angles}


                    coordinated_metals[metal_res] = coordinated_atoms_list#, 'angles': angles}
    coordinated_atoms_ids = {}
    if coordinated_metals:
        for metal, atoms_list in coordinated_metals.iteritems():
            metal_id = "{} {} {}".format(metal.getResname(), metal.getChid(), metal.getResnum())
            atoms_ids = [["{} {} {} {}".format(at.getResnum(), at.getResname(), at.getChid(), at.getName()),
                          calcDistance(metal, at)[0]] for at in atoms_list]
            if len(atoms_list) in [x[1] for x in coordination_geometries.itervalues()]:
                coordinated_atoms_ids[metal_id] = atoms_ids
            print "     * The metal atom {0} has the following atoms within coordination " \
                  "distance:\n{1}".format(metal_id, "\n".join(['       * {0}'.format(x[0]) for x in atoms_ids]))
            angles = [[at, metal, at2, calcAngle(at, metal, at2)[0]]
                      for idx,at in enumerate(atoms_list) for at2 in atoms_list[idx + 1:]]
            if len(atoms_list) <= 4:
                print "      * Checking for a tetrahedric coordination for the atom."
                found_conformation = CheckConformation(angles, 'tetrahedric')
                if not found_conformation:
                    print "        * WARNING: The angles are too distorted to ascertain this configuration. CHECK IT manually"
            elif 4 < len(atoms_list) <= 6: #or not found_conformation:
                print "      * WARNING: Checking for an octahedric coordination for the atom."
                found_conformation = CheckConformation(angles, 'octahedric')
                if not found_conformation:
                    print "        * The angles are too distorted to ascertain this configuration. CHECK IT manually"
            else:
                print "      * The metal doesn't have a coordination we can validate."


    else:
        print "  * There are no coordinated metals."

    return coordinated_metals


def CheckStructure(initial_structure, gaps={}, no_gaps={}, charge_terminals=False, remove_missing_ter=False,
                   debug=False):
    residues2fix = {}
    crosslinked_cysteines, charged_cysteines = CheckCysteines(initial_structure)
    residues2remove = {}
    metals2coordinate = CheckMetalsCoordination(initial_structure)
    residues_without_template = []
    for chain in initial_structure.iterChains():
        if chain.getChid() in gaps.keys():
            gaps_e = [x[0] for x in gaps[chain.getChid()]]
            gaps_b = [x[1] for x in gaps[chain.getChid()]]
        else:
            gaps_e = []
            gaps_b = []
        initial_residue, final_residue = FindInitialAndFinalResidues(chain)
        residues2ignore = []
        for residue in chain.iterResidues():
            resname = residue.getResname().strip()
            resnum = residue.getResnum()
            residue_atomnames = list(residue.getNames())
            if resname in supported_metals:
                continue
            if resname in supported_aminoacids:
                if charge_terminals:
                    if resnum == initial_residue or resnum in gaps_b:
                        zmatrix = ZMATRIX(resname + 'B')
                    elif resnum == final_residue or resnum in gaps_e:
                        zmatrix = ZMATRIX(resname + "E")
                    else:
                        zmatrix = ZMATRIX(resname)
                else:
                    if resname == "NMA":
                        zmatrix = ZMATRIX(resname + "E")
                    elif resname == "ACE":
                        zmatrix = ZMATRIX(resname + "B")
                    else:
                        zmatrix = ZMATRIX(resname)
            else:
                try:
                    zmatrix = ZMATRIX(resname)
                except:
                    print "UP {}".format(resname)
            if zmatrix.Name is None:
                print "  * The residue {} {} doesn't have a template, so it won't be checked.".format(resname, resnum)
                residues_without_template.append([resname, resnum, chain.getChid()])
                continue
            residue_atomnames.sort()
            if sorted(zmatrix.AtomNames) != residue_atomnames:
                atoms2delete = []
                atoms2add = []
                atoms2modify = []
                missing_atoms = [atom_name for atom_name in zmatrix.AtomNames if atom_name not in residue_atomnames]
                over_atoms = [atom_name for atom_name in residue_atomnames if atom_name not in zmatrix.AtomNames]
                if over_atoms:
                    if charge_terminals:
                        if resnum == final_residue or resnum in gaps_e:
                            if "HXT" in over_atoms:
                                over_atoms.pop(over_atoms.index('HXT'))
                                if 'OXT' in missing_atoms:
                                    atoms2modify.append(['HXT', 'OXT'])  # This changes the atom name
                                    missing_atoms.pop(missing_atoms.index('OXT'))
                                else:
                                    atoms2delete.append('HXT')
                        elif resnum == initial_residue or resnum in gaps_b:
                            if "H1" in over_atoms:
                                over_atoms.pop(over_atoms.index('H1'))
                            if "H2" in over_atoms:
                                over_atoms.pop(over_atoms.index('H2'))
                        if over_atoms:
                            print "   The residue {} {} has more atoms than the zmatrix." \
                                  " The extra atoms are:{}\n PELE won't work, review " \
                                  "them".format(resname, resnum, ",".join(over_atoms))
                    else:
                        if resnum == initial_residue or resnum in gaps_b:
                            if "H1" in over_atoms:
                                over_atoms.pop(over_atoms.index('H1'))
                                if 'H' in missing_atoms:
                                    missing_atoms.pop(missing_atoms.index('H'))
                                    atoms2modify.append(['H1', 'H'])  # This changes the atom name
                                else:
                                    atoms2delete.append('H1')
                            if "H2" in over_atoms:
                                over_atoms.pop(over_atoms.index('H2'))
                                atoms2delete.append('H2')
                            if "H3" in over_atoms:
                                over_atoms.pop(over_atoms.index('H3'))
                                atoms2delete.append('H3')
                        elif resnum == final_residue or resnum in gaps_e:
                            if "HXT" in over_atoms:
                                over_atoms.pop(over_atoms.index('HXT'))
                                atoms2delete.append('HXT')
                            if "OXT"in over_atoms:
                                over_atoms.pop(over_atoms.index('OXT'))
                                atoms2delete.append('OXT')
                        if over_atoms:
                            print "   The residue {} {} has more atoms than the zmatrix. The extra atoms are:{}\n " \
                                  "PELE won't work, review them".format(resname, resnum, ",".join(over_atoms))
                    # residues2fix["{} {} {}".format(resname, resnum, residue.getChid())] = {'delete': atoms2delete}
                    if debug:
                        print 'over_atoms: ', over_atoms
                        print "ZMATRIX: {}\n zmat atoms:{}".format(zmatrix.Name, sorted(zmatrix.AtomNames))
                        print " resi atoms:", sorted(residue_atomnames)
                if missing_atoms:
                    if "CA" in missing_atoms or "N" in missing_atoms or "O" in missing_atoms or "C" in missing_atoms:
                        print "  The residue {} {} is missing one or more heavy atoms in the backbone. " \
                              "It cannot be fixed.".format(resname, resnum, residue.getChid())
                        if remove_missing_ter:
                            if resnum in [initial_residue, final_residue]:
                                print " INFO: In structure {} the residue {} {} in chain {} will be eliminated".format(
                                    initial_structure.getTitle(), residue.getResname(),
                                    resnum, residue.getChid())
                                try:
                                    residues2remove[residue.getChid()]
                                except KeyError:
                                    residues2remove[residue.getChid()] = [resnum]
                                else:
                                    residues2remove[residue.getChid()].append(resnum)
                                if resnum == initial_residue:
                                    next_res_resnum = resnum + 1
                                    next_residue = chain.getResidue(next_res_resnum)
                                    next_res_resname = next_residue.getResname()
                                    if charge_terminals:
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
                                        key = "{} {} {}".format(next_res_resname, next_res_resnum, residue.getChid())
                                        residues2fix[key]['add'] = atoms2add
                                        residues2fix[key]['delete'] = []
                                        residues2fix[key]['modify'] = []
                                elif resnum == final_residue:
                                    prev_res_resnum = resnum - 1
                                    prev_residue = chain.getResidue(prev_res_resnum)
                                    if prev_residue is None:
                                        print "There's no residue with number {} in chain {}".format(prev_res_resnum,
                                                                                                     chain.getChid())
                                    elif charge_terminals:
                                        prev_res_resname = prev_residue.getResnames()[0]
                                        atoms2add = {"OXT"}
                                        key = "{} {} {}".format(prev_res_resname, prev_res_resnum, residue.getChid())
                                        try:
                                            residues2fix[key]
                                        except KeyError:
                                            residues2fix[key]['add'] = atoms2add
                                            residues2fix[key]['delete'] = []
                                            residues2fix[key]['modify'] = []
                                        else:
                                            residues2fix[key]['add'] += atoms2add
                        else:
                            print "  The program will be interrupted.".format(resname, resnum)
                            print "     This are the missing atoms: {}".format(missing_atoms)
                            sys.exit()
                    else:
                        if "H" in missing_atoms and resnum in gaps_b and resnum not in residues2ignore:
                            maestro_terminal_H = ["H1", "H2", "H11", "H22"]
                            possible_atoms2change = [atom_name for atom_name in maestro_terminal_H
                                               if atom_name in residue.getNames()]
                            if possible_atoms2change:
                                atomname_to_use = [0]
                                if atomname_to_use:
                                    atom = residue.getAtom(atomname_to_use)
                                    atom.setName('H')
                                    missing_atoms.pop(missing_atoms.index('H'))
                        if resname == "CYS":
                            key = "{}_{}".format(residue.getChid(), resnum)
                            if key in crosslinked_cysteines + charged_cysteines:
                                if key in charged_cysteines:
                                    residue.setResname('CYT')
                                elif "HG" in missing_atoms:
                                    missing_atoms.pop(missing_atoms.index('HG'))
                        for atom_name in missing_atoms:
                            if atom_name[0] != "H" and atom_name[0].isalpha() and atom_name != 'OXT':
                                print "  The residue {} {} is missing the heavy atom {}.\n" \
                                      "  All the atoms depending on this atom will" \
                                      " be placed according to the zmatrix.".format(resname, resnum, atom_name)
                                atoms2add = zmatrix.GetAllChildrenAtoms(atom_name)
                                atoms2add = set(atoms2add).union(set(missing_atoms))
                            elif resnum not in residues2ignore:
                                atoms2add = set(missing_atoms)
                key = "{} {} {}".format(resname, resnum, residue.getChid())
                try:
                    residues2fix[key]
                except KeyError:
                    if atoms2add or atoms2delete or atoms2modify:
                        residues2fix[key] = {'add': atoms2add, 'delete': atoms2delete, 'modify': atoms2modify}
                else:
                    residues2fix[key]['add'] += atoms2add,
                    residues2fix[key]['delete'] += atoms2delete
                    residues2fix[key]['modify'] += atoms2modify
    # print residues2fix, residues2remove
    return residues2fix, residues2remove, metals2coordinate, residues_without_template


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
                # if residue.getResnum() > previous_residue_number + 1:
                if residue.getResname() in supported_aminoacids:
                    if residue.getResname() == 'ACE':
                        try:
                            gaps[chain_id]
                        except KeyError:
                            gaps[chain_id] = []
                        gaps[chain_id].append([previous_residue_number, residue.getResnum()])
                    else:
                        current_residue_N = residue.getAtom("N")
                        previous_residue = chain.getResidue(previous_residue_number)
                        previous_residue_C = previous_residue.getAtom('C')
                        if current_residue_N is not None and previous_residue_C is not None:
                            distance = calcDistance(current_residue_N, previous_residue_C)
                            if distance < 2.5:
                                try:
                                    not_gaps[chain_id]
                                except KeyError:
                                    not_gaps[chain_id] = []
                                not_gaps[chain_id].append([previous_residue_number, residue.getResnum()])
                            else:
                                try:
                                    gaps[chain_id]
                                except KeyError:
                                    gaps[chain_id] = []
                                gaps[chain_id].append([previous_residue_number, residue.getResnum()])
                        elif current_residue_N is None:
                            print "   * There's a problem with residue {} {} {} it" \
                                  " doesn't have the N atom".format(residue.getResname(),
                                                                    residue.getResnum(),
                                                                    residue.getChid())
                        elif previous_residue_C is None:
                            print "   * There's a problem with residue {} {} {} it doesn't have the C atom".format(
                                previous_residue.getResname(),
                                previous_residue.getResnum(),
                                previous_residue.getChid())
                if residue.getResnum() == previous_residue_number:
                    if residue.getIcode() == '':
                        if residue.hetero is not None:
                            residue.setResnum(previous_residue_number + 1)
                        else:
                            print '   * WARNING: There are two residues with the same resnum in the same chain.' \
                                  'The residue {} {} {} '.format(residue.getChid(),
                                                                 residue.getResnum(),
                                                                 residue.getResname())
                elif residue.getResnum() < previous_residue_number + 1:
                    if residue.hetero is not None:
                        residue.setResnum(previous_residue_number + 1)
                    else:
                        # print residue.getResnum(), previous_residue_number
                        print "   * WARNING! The next residue has a lower resnum than the previous one in the chain??"
                        print "    * CHAIN: {} RESIDUE: {}  {} {}".format(residue.getChid(),
                                                                    residue.getResnum(),
                                                                    residue.getChid(),
                                                                    residue.getResname())
                        print "    * previous: RESIDUE {} {}  {}".format(previous_residue.getResname(),
                                                                   previous_residue.getResnum(),
                                                                   previous_residue.getChid())
                        sys.exit('This should never happen, please review your structure.')
                else:
                    pass
            previous_residue_number = residue.getResnum()
    return gaps, not_gaps
