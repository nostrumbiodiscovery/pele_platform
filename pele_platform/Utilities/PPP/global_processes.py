import sys
from argparse import ArgumentParser, RawTextHelpFormatter
from re import match
from pele_platform.Utilities.PPP.global_variables import supported_aminoacids, input_keywords, aminoacids_1letter, aminoacids_3letter
from pele_platform.Utilities.PPP import parameters_help

from prody import parsePDB

__author__ = 'jelisa'


def ParseMutations(unprocessed_mutation, structure_filename):
    """"""
    # TODO: This fuction should be able to allow the user to specify
    # regions and selections in the input instead of just one residue.
    processed_mutations = []
    chain = ""
    for mutation in unprocessed_mutation:
        splitted_mutation = mutation.split()
        if "residue" not in mutation:
            if len(splitted_mutation) == 4:
                original_resname, resnum, chain, final_resname = [value.upper() for value in splitted_mutation]
            elif len(splitted_mutation) == 3:
                original_resname, resnum, final_resname = [value.upper() for value in splitted_mutation]
            else:
                selection_string_raw, final_resname_raw = mutation.split("TO")
                structure = parsePDB(structure_filename)
                #TODO: Talk with fatima about the selections they'll be doing.
                if "within" in selection_string_raw:
                    selection_center = selection_string_raw.split("of")[-1]
                    selection_string = "({}) and (not {}) and (not resname HOH)".format(selection_string_raw,
                                                                                        selection_center)
                    mutations2add = []
                    selection = structure.select(selection_string).copy()
                    if selection is None:
                        print("Something is wrong with the selection string it isn't selecting anything. Review it!")
                        print(" Selection string: {}".format(selection_string))
                        print("The program will be interrupted.")
                        sys.exit()
                    for residue in selection.iterResidues():
                        original_resname = residue.getResname()
                        resnum = residue.getResnum()
                        chain = residue.getChid()
                        final_resname = final_resname_raw
                        # mutations2add.append(mutation)
                else:
                    selection = selection_string_raw.split()
                    original_resname = selection[selection.index("resname") + 1].upper()
                    resnum = selection[selection.index("resnum") + 1]
                    if "chain" in selection:
                        chain = selection[selection.index("chain") + 1].upper()
                    else:
                        chain = ""
                    final_resname = final_resname_raw.strip().upper()
        else:
            if len(splitted_mutation) != 4:
                try:
                    int(splitted_mutation[splitted_mutation.index("residue") + 2])
                except ValueError:
                    try:
                        int(splitted_mutation[splitted_mutation.index("residue") + 1])
                    except ValueError:
                        print("Your mutation string is wrong." \
                              " It should look like: 'residue XXX N to YYY' with N being the resnum and XXX, YYY being supported aminoacids")
                        sys.exit()
                    else:
                        original_resname = splitted_mutation[splitted_mutation.index("residue") + 2].upper()
                        resnum = splitted_mutation[splitted_mutation.index("residue") + 1]
                else:
                    original_resname = splitted_mutation[splitted_mutation.index("residue") + 1].upper()
                    resnum = splitted_mutation[splitted_mutation.index("residue") + 2]
                try:
                    chain = splitted_mutation[splitted_mutation.index("chain") + 1]
                except ValueError:
                    if splitted_mutation[splitted_mutation.index("residue") + 3] != "to":
                        chain = splitted_mutation[splitted_mutation.index("residue") + 3]
                    else:
                        chain = ""
                final_resname = splitted_mutation[splitted_mutation.index('to') + 1].upper()
        mutation = {'ini_resname': original_resname, "resnum": resnum,
                    "chain": chain, "fin_resname": final_resname}
        if mutation["fin_resname"] in supported_aminoacids:
            processed_mutations.append(mutation)
            print('a')
        elif mutation["fin_resname"].lower() in input_keywords.keys():
            keyword = mutation["fin_resname"].lower()
            histidines_equivalents = ["HID", "HIE", "HIP", "HIS"]
            if original_resname in histidines_equivalents:
                exclude_histidines = True
            else:
                exclude_histidines = False
            for aminoacid in input_keywords[keyword]:
                if exclude_histidines and aminoacid in histidines_equivalents:
                    continue
                elif aminoacid == original_resname:
                    continue
                else:
                    mutation = {'ini_resname': original_resname, "resnum": resnum,
                                "chain": chain, "fin_resname": aminoacid}
                    processed_mutations.append(mutation)
        else:
            print('c', mutation["fin_resname"], input_keywords.keys())
    return processed_mutations


def ParseArguments():
    parser = ArgumentParser(description=parameters_help.program_description,
                            formatter_class=RawTextHelpFormatter)
    parser.add_argument('-input_pdb', '-ipdb', required=True, nargs='*', help=parameters_help.ipdb_description)
    parser.add_argument('-output_pdb', '-opdb', default=[], nargs="+", help=parameters_help.opdb_description)
    parser.add_argument('-mutation', '-mut', type=str, default='', nargs='+',
                        help=parameters_help.mutation_description)
    parser.add_argument("-mutants_from_file", '-mut_file', default=False,
                        help=parameters_help.mutants_from_file_description)
    parser.add_argument("-mutant_multiple", action="store_true",
                        help=parameters_help.mutant_multiple_description)
    parser.add_argument("-charge_terminals", action="store_true",
                        help=parameters_help.charge_terminals_description)
    parser.add_argument("-no_gaps_ter", action="store_true", default=False,
                        help=parameters_help.no_gaps_ter_description)
    parser.add_argument("-pdb_resolution", default=2.5, type=float,
                        help=parameters_help.pdb_resolution_description)
    parser.add_argument("-make_unique", default=False,
                        help=parameters_help.make_unique_description)
    parser.add_argument("-remove_terminal_missing", default=False, action="store_true",
                        help=parameters_help.remove_terminal_missing_description)
    # parser.add_argument('-force_field', '-ff', type=str, default="OPLS2005",
    #                     choices=["OPLS2005", "OPLS2001", "AMBER99sbBSC0"],
    #                     help=parameters_help.ff_description)
    args = parser.parse_args()

    if len(args.input_pdb) > 1 and (args.mutation != '' or args.mutants_from_file):
        print("This options are incompatible, when you want to do mutations you should use only one structure at time.")
        return None
    if len(args.input_pdb) == 1:
        args.input_pdb = args.input_pdb[0]
    if args.mutation:
        args.mutation = ParseMutations(args.mutation, args.input_pdb)
    if args.mutants_from_file:
        filein = open(args.mutants_from_file, 'r')
        mutations_text = filein.readlines()
        filein.close()
        args.mutation = []
        args.output_pdb = []
        for line in mutations_text:
            line = line.strip()
            if line == '':
                continue
            output_file, mutation = line.split('\t')
            splitted_mutation = mutation.split()
            if len(splitted_mutation) == 4:
                args.mutation.append(
                    {'ini_resname': splitted_mutation[0].upper(), "resnum": splitted_mutation[1].upper(),
                     "chain": splitted_mutation[2].upper(), "fin_resname": splitted_mutation[3].upper()})
            elif len(splitted_mutation) == 3:
                args.mutation.append(
                    {'ini_resname': splitted_mutation[0].upper(), "resnum": splitted_mutation[1].upper(),
                     "fin_resname": splitted_mutation[3].upper()})
    output_files_number = len(args.output_pdb)
    if output_files_number == 0:
        args.output_pdb = []
        args.output_pdb.append(args.input_pdb[:-4] + "_processed" + args.input_pdb[-4:])
    elif output_files_number == 1:
        if len(args.mutation) > 1:
            general_name = args.mutation[0].split('.pdb')[0]
            for mutation in args.mutation:
                mutation_id = aminoacids_1letter[aminoacids_3letter.index(mutation["ini_resname"])]
                if mutation["chain"]:
                    mutation_id += mutation["chain"]
                mutation_id += mutation["resnum"]
                mutation_id += aminoacids_1letter[aminoacids_3letter.index(mutation["fin_resname"])]
                output_name = "{}_{}.pdb".format(general_name, mutation_id)
                args.output_pdb.append(output_name)


    return args


def FindInitialAndFinalResidues(structure):
    """
    This function finds the first and the last residue of a protein.
    """

    initial_residue, final_residue = 0, 0
    first_aa = False
    for res in structure.iterResidues():
        resnum = res.getResnum()
        resname = res.getResname()
        if resname in supported_aminoacids:
            if not first_aa:
                initial_residue = resnum
                first_aa = True
            elif resnum > final_residue:
                final_residue = resnum

    return initial_residue, final_residue


def PDBwriter(output_file_name, structure, make_unique, residues2remove, no_ter_mkar_for_gaps, no_proteic_ligand=None,
              gaps={}, no_gaps={}, mid_chain_nonstd_residue=[]):
    supported_aminoacids.extend(mid_chain_nonstd_residue)
    if '.pdb' not in output_file_name:
        output_file_name += '.pdb'
    outfile = open(output_file_name, 'w')
    alt_loc = " "  # It should always be an empty string.
    insertion_code = " "  # It should always be an empty string.
    occupancy = 1.0  # It should be always 1
    serial = 1
    for chain in structure.iterChains():
        ter = False
        resnums = [residue.getResnum() for residue in chain]
        if chain.getChid() in gaps.keys():
            gaps_e = [x[0] for x in gaps[chain.getChid()]]
        else:
            gaps_e = []
        if chain.getChid() == make_unique and no_proteic_ligand is not None:
            change_names = True
        else:
            change_names = False
        ligand_possible_atoms = {"C": 1, "N": 1, "O": 1, "H": 1, "F": 1,
                                 "S": 1, "P": 1, "BR": 1, "I": 1, "CL": 1,
                                 "FE": 1, "CA": 1}

        for index, residue in enumerate(chain.iterResidues()):
            ter = False
            resnum = residue.getResnum()
            eliminate_hxt = False
            eliminate_h = False
            if residue.getChid() in residues2remove.keys():
                if resnum in residues2remove[residue.getChid()]:
                    continue
            if resnum in gaps_e:
                ter = True
                if no_ter_mkar_for_gaps:
                    ter = False
                if " HXT" in residue.getNames():
                    eliminate_hxt = True
                elif " H1 " in residue.getNames() or " H2 " in residue.getNames():
                    eliminate_h = True
            chain_id = residue.getChid()
            resname = residue.getResname()
            if resname == "NMA":
                ter = True
                atom_hetatm = "ATOM"
            elif resname in supported_aminoacids:
                atom_hetatm = "ATOM"
                try:
                    next_residue = chain.select("resnum `{}`".format(resnums[index+1]))
                except IndexError:
                    pass
                else:
                    if next_residue.getResnames()[0] not in supported_aminoacids:
                        ter = True
            else:
                atom_hetatm = "HETATM"
                ter = True
            for atom in residue.iterAtoms():
                atom_name = atom.getName()
                if atom_name == ' HXT' and eliminate_hxt:
                    continue
                elif atom_name in [' H1 ', ' H2 '] and eliminate_h:
                    continue
                if change_names:
                    if "Ca" in atom_name:
                        raw_atom_name = atom_name.strip()[:2]
                    elif "CA" in atom_name:
                        if atom_name[:2] == "CA":
                            raw_atom_name = atom_name.strip()[:2]
                        else:
                            raw_atom_name = atom_name.strip()[0]
                    else:
                        raw_atom_name = atom_name.strip()
                        if len(raw_atom_name) > 1:
                            # print "'{}'".format(atom_name)
                            raw_atom_name = raw_atom_name.upper()
                            if raw_atom_name[:2] in ligand_possible_atoms.keys():
                                raw_atom_name = raw_atom_name[:2]
                            elif raw_atom_name not in ligand_possible_atoms.keys():
                                raw_atom_name = raw_atom_name[0]
                                if raw_atom_name not in ligand_possible_atoms.keys():
                                    print("INVALID ATOM NAME in the ligand file: {}".format(raw_atom_name))
                        elif raw_atom_name not in ligand_possible_atoms.keys():
                            print("INVALID ATOM NAME in the ligand file: {}".format(raw_atom_name))
                    atom_name = raw_atom_name + str(ligand_possible_atoms[raw_atom_name])
                    if len(atom_name) == 1:
                        formated_atom_name = " " + atom_name + "  "
                    elif len(atom_name) == 2:
                        formated_atom_name = " " + atom_name + " "
                    elif len(atom_name) == 3:
                        formated_atom_name = " " + atom_name
                    ligand_possible_atoms[raw_atom_name] += 1
                else:
                    formated_atom_name = atom_name
                x, y, z = atom.getCoords()
                element = atom.getElement()
                segment = atom.getSegname()
                b_factor = atom.getBeta()
                if atom.getCharge() is None:
                    charge = " "
                else:
                    charge = atom.getCharge()
                pdbstring = "{1:6}{2:5}{0:1}{3:4}{4:1}{5:3}{0:1}{6:1}{7:4}{8:1}{0:3}" \
                            "{9:8.3f}{10:8.3f}{11:8.3f}{12:6.2f}{13:6.2f}{0:6}{14:<4}{15:>2}{16:2}\n"
                pdbstring = pdbstring.format('', atom_hetatm, serial, formated_atom_name, alt_loc,
                                             resname, chain_id, resnum, insertion_code, x, y,
                                             z, occupancy, b_factor, segment,
                                             element, charge)
                outfile.write(pdbstring)
                serial += 1
            if ter:
                outfile.write("TER\n")
        if not ter:
            outfile.write("TER\n")
    outfile.close()
    return


def RenumberStructure(initial_structure, debug=False):
    renumbered_structure = initial_structure.copy()
    for chain in renumbered_structure.iterChains():
        residue_number = 1
        chain_id = chain.getChid()
        for residue in chain.iterResidues():
            if residue.getResnum() == debug:
                print(residue.getResnum(), residue.getResname(), residue_number)
            residue.setResnum(residue_number)
            residue.setIcode('')
            residue_number += 1
    return renumbered_structure
