import copy

from prody import calcDistance

from checks_module import CheckClashes
from coordinates_module import ChangeResidueCoordinates
from global_processes import FindInitialAndFinalResidues
from global_variables import protein_atomnames_dictionary, supported_aminoacids, supported_metals
from program_own_classes import ROTAMERLIB

__author__ = 'jelisa'


def SolveClashes(initial_structure, initial_clashes, mutation, zmatrix, initial_residue_number, final_residue_number):
    print "Trying to solve the clashes."
    dihedral2check = 0
    initial_part_of_the_protein = initial_structure.select("resnum < {}".format(mutation["resnum"])).copy()
    final_part_of_the_protein = initial_structure.select("resnum > {}".format(mutation["resnum"])).copy()
    initial_residue = initial_structure.select("resnum {}".format(mutation["resnum"])).copy()
    resname = initial_residue.getResnames()[0]

    if resname in ["HIE", "HID", "HIP"]:
        resname = "HIS"
    elif resname in ["LYN"]:
        resname = "LYS"
    residue_rotamer_lib = ROTAMERLIB(resname)
    maximum_number_of_trials = residue_rotamer_lib.NumberOfDihedrals
    initial_number_of_clashes = 0
    for values in initial_clashes.values():
        initial_number_of_clashes += len(values)
    print "Starting number of heavy atoms clashes: {}".format(len(initial_clashes.keys()))
    print "Starting number of total clashes: {}".format(initial_number_of_clashes)

    zmatrix.ComputeDeltaFi()
    while maximum_number_of_trials > 0:
        maximum_number_of_trials -= 1
        new_residue = ChangeResidueCoordinates(initial_residue.copy(), zmatrix, residue_rotamer_lib, dihedral2check)
        if mutation["fin_resname"] == "PRO":
            CD_N_distance = calcDistance(new_residue.select("name N"), new_residue.select("name CD"))
            if CD_N_distance > 1.5:
                dihedral2check += 1
                continue
        dihedral2check += 1
        mutated_structure = initial_part_of_the_protein + new_residue + final_part_of_the_protein
        new_clashes = CheckClashes(mutated_structure, mutation, zmatrix, initial_residue_number, final_residue_number)
        current_number_of_clashes = 0
        for values in new_clashes.values():
            current_number_of_clashes += len(values)
        if new_clashes == {}:
            initial_residue = new_residue
            break
        elif len(new_clashes.keys()) < len(initial_clashes.keys()):

            print 'The number of clashing heavy atoms in the residue now is: {} and clashes with {} atoms'.format(
                len(new_clashes.keys()),
                current_number_of_clashes)
            initial_residue = new_residue.copy()
            initial_number_of_clashes = current_number_of_clashes
            initial_clashes = new_clashes
        elif len(new_clashes.keys()) == len(initial_clashes.keys()):
            if current_number_of_clashes < initial_number_of_clashes:
                print "The number of clashing heavy atoms is still the same but now they clash with {} atoms.".format(
                    current_number_of_clashes)
                initial_residue = new_residue.copy()
                initial_number_of_clashes = current_number_of_clashes
    return initial_part_of_the_protein + initial_residue + final_part_of_the_protein


def FixAtomNames(initial_structure, gaps={}, no_gaps={}, debug=False):

    correct_structure = initial_structure.copy()

    for chain in correct_structure.iterChains():
        # for chain in structure_without_waters.iterChains():
        if chain.getChid() in gaps.keys():
            gaps_residues = [y for x in gaps[chain.getChid()] for y in x]
            gaps_e = [x[0] for x in gaps[chain.getChid()]]
            gaps_b = [x[1] for x in gaps[chain.getChid()]]
        else:
            gaps_residues = []
            gaps_b = []
            gaps_e = []
        if chain.getChid() in no_gaps.keys():
            no_gaps_residues = [y for x in no_gaps[chain.getChid()] for y in x]
        else:
            no_gaps_residues = []
        initial_res, final_res = FindInitialAndFinalResidues(chain)
        if debug:
            print "working with chain {}".format(chain.getChid())
        for residue in chain.iterResidues():
            tmp_dictio = copy.deepcopy(protein_atomnames_dictionary)
            resname = residue.getResname()
            try:
                possible_names = tmp_dictio[resname]
            except KeyError:
                print '   * The residue {} is not an aa nor a water.'.format(resname)
                possible_names = [["CL"], ["CU"], ["FE1"], ["FE2"], ["ZN"], ["MG"]]
                heteroatom = True
            else:
                if resname == 'HOH':
                    heteroatom = True
                else:
                    heteroatom = False
            if residue.getResnum() == initial_res or residue.getResnum() in gaps_b:
                # This magic number 0 comes from the fact that for all the aminoacids
                # except the proline the first set of atoms in the protein_atomnames_dictionary
                # should have the hydrogen H from the nitrogen, but if the residue is the
                # first one in the protein it should be changed into H3 which belongs to
                # the "END" keyword in the names dictionary.
                if debug:
                    print 'the initial residue is: {} and the residue number: {}'.format(initial_res,
                                                                                         residue.getResnum())
                if 'H' in possible_names[0]:
                    possible_names.pop(0)
            if residue.getResnum() in [initial_res, final_res] + gaps_residues:
                # The other in this step is really important it should always be the
                # possible and then the ending possibilities, otherwise it will create
                # problems for the ending residues.
                possible_names = possible_names + tmp_dictio["END"]

            # if debug and residue.getResnum() == debug[0] and chain.getChid() == debug[1]:
            #     print residue.getResname(), residue.getResnum(), residue.getChid(), residue.getNames()
            for atom in residue.iterAtoms():
                atom_found = False
                atom_name = atom.getName()
                # if debug and residue.getResnum() == debug[0] and chain.getChid() == debug[1]:
                #     print 'pos', atom_name, possible_names
                for atoms in possible_names:
                    if atom_name in atoms:
                        old_atom_name = atom_name
                        atom_name = atoms[0]
                        atom_found = True
                        possible_names.pop(possible_names.index(atoms))
                        if debug and atom.getResnum() == debug[0] and chain.getChid() == debug[
                            1]:  # .split()[0] and atom.getResnum() == debug.split()[1]:
                            print 'a', residue.getResnum(), old_atom_name, atom_name, final_res
                            # print 'a', possible_names

                        break
                if not atom_found and not heteroatom:
                    if atom_name in ['HXT', 'H1', "H2"]:
                        if residue.getResnum() in gaps_residues or residue.getResnum() in [initial_res, final_res]:
                            pass
                        elif residue.getResnum() in no_gaps_residues:
                            print "   * The residue {} won't be considered as a gap, if it really is one," \
                                  " let the developer know".format("{2} {1} {0}".format(residue.getResnum(),
                                                                                        residue.getChid(),
                                                                                        residue.getResname()))
                    else:
                        print "   * The Atom {} from residue {} {} {} doesn't have a valid name.".format(atom_name,
                                                                                                         resname,
                                                                                                         atom.getChid(),
                                                                                                         atom.getResnum())
                atom.setName(atom_name)
                # if debug: break
    return correct_structure


def WritingAtomNames(initial_structure):
    correct_names_structure = initial_structure.copy()
    for atom in correct_names_structure.iterAtoms():
        atom_name = atom.getName()
        if atom.getResname() in supported_metals:
            atom_name = "{0}  ".format(atom_name)
        elif atom_name == "1HW":
            atom_name = "1HW "
        elif atom_name == "2HW":
            atom_name = "2HW "
        elif atom_name == "FE1":
            atom_name = "FE1 "
        elif atom_name == "FE2":
            atom_name = "FE2 "
        elif atom.getResname() == "NA":
            atom_name = "NA  "
        elif atom_name == "CL":
            atom_name = "CL  ".format(atom_name)
        elif len(atom_name) == 1:
            atom_name = " " + atom_name + "  "
        elif len(atom_name) == 2:
            atom_name = " " + atom_name + " "
        elif len(atom_name) == 3:
            atom_name = " " + atom_name
        atom.setName(atom_name)
    return correct_names_structure


def FixStructureResnames(initial_structure, ligand_chain=False):
    """
    This function changes the residues names so they match with the names required by PELE.
    """
    structure = initial_structure.copy()
    for res in structure.iterResidues():
        resname = res.getResname()
        ligand_structure = None
        # print ligand_chain == res.getChid(), resname in supported_aminoacids
        if ligand_chain == res.getChid() and resname in supported_aminoacids:
            ligand_structure = structure.chain_Z.copy()
        if ligand_structure is not None and ligand_structure.hetero is not None:
            print "INFO: Renaming the ligand to LIG in structure {}".format(initial_structure.getTitle())
            resname = "LIG"
            res.setResnum(1)
        if resname == "HIS":
            if res.numAtoms() in [10, 16]:
                resname = "HIE"
            elif res.numAtoms() == 19:
                if ("HE2" in res.getNames() and "HD1" in res.getNames()) or\
                        ("HNE" in res.getNames() and "HND" in res.getNames()):
                    resname = "HIP"
                elif "HE2" in res.getNames() or "HNE" in res.getNames():
                    resname = "HIE"
                elif "HD1" in res.getNames() or "HND" in res.getNames():
                    resname = "HID"
            elif res.numAtoms() == 18:
                if "HE2" in res.getNames() and "HD1" in res.getNames() or\
                        ("HNE" in res.getNames() and "HND" in res.getNames()):
                    resname = "HIP"
                elif "HE2" in res.getNames() or "HNE" in res.getNames():
                    resname = "HIE"
                elif "HD1" in res.getNames() or "HND" in res.getNames():
                    resname = "HID"
                else:
                    print "Make sure the histidine is correctly protonated and then, only then, explain the bug" \
                          " to the the developer."
            elif res.numAtoms() == 17:
                if "HE2" in res.getNames() or "HNE" in res.getNames():
                    resname = "HIE"
                elif "HD1" in res.getNames() or "HND" in res.getNames():
                    resname = "HID"
            else:
                print "  * WARNING: The hystidine {} has an incorrect number " \
                      "of atoms. PELE won't work correctly".format("{} {}".format(res.getResnum(), res.getChid(),))
        elif resname == "LYS":
            if res.numAtoms() == 21:
                resname = "LYN"
            elif "OXT" in res.getNames() and res.numAtoms() == 22:
                resname = "LYN"
            elif "H1" in res.getNames() and res.numAtoms() == 23 and "HZ3" not in res.getNames():
                    resname = "LYN"
        elif resname == "GLU":
            if res.hydrogen is not None:
                if res.numAtoms() == 16 and res.hydrogen.numAtoms() == 7 and "H1" not in res.getNames() and \
                                "HXT" not in res.getNames():
                    resname = "GLH"
                elif res.numAtoms() == 18 and res.hydrogen.numAtoms() == 9:
                    resname = "GLH"
                elif res.numAtoms() == 17 and res.oxygen.numAtoms() == 4:
                    resname = "GLH"
        elif resname == "ASP":
            if res.numAtoms() == 13 and "OXT" not in res.getNames() and "HXT" not in res.getNames():
                resname = "ASH"
            elif "OXT" in res.getNames() and res.numAtoms() == 14:
                resname = "ASH"
            elif "H1" in res.getNames() and res.numAtoms() == 15:
                resname = "ASH"
        elif resname == "CL":
            resname = " CL"
        res.setResname(resname)
    return structure


def MoveBackbone(initial_structure, zmatrix, initial_part_of_the_protein, final_part_of_the_protein):
    return
