__author__ = 'jelisa'

program_description = "This program will substitute plop and ppp for PELE. It'll check the structure\n" \
                      " the atom names in the future, for now it does mutations and checks-formats the\n" \
                      " names of the residues, produces a .pdb file that PELE can read without problems."
ipdb_description = "This parameter is mandatory and specifies the path to\n" \
                   " the .pdb file with the structure to mutate.\n" \
                   "This structure should have all the atoms present,\n" \
                   " including the hydrogens."
opdb_description = "This specifies the path to write the final pdb.\n" \
                   "If it isn't specified the program will generate\n" \
                   " a pdb with the same name as the initial pdb but\n" \
                   " with the suffix _processed."
mutation_description = "This parameter specifies the desired mutation(s).\n" \
                       "It should be a string wiht the desired mutation(s)\n" \
                       " in the following format:\n" \
                       "  'residue XXX N to YYY'\n" \
                       "  where:\n" \
                       "    - XXX: the original residue\n" \
                       "    - N: the number of the residue\n" \
                       "    - YYY: the desired aminoacid.\n" \
                       "The names of the residues should be the three letters\n" \
                       " code for the aminoacid.\n" \
                       "EXCEPTIONS: the hystidine you should use one of the \n" \
                       " following names: HID, HIP, HIE"
ff_description = "This specifies the force-field to use. The supported force-fields are:\n" \
                 " - OPLS2005 (default)\n" \
                 " - OPLS2001\n" \
                 " - AMBER99sbBSC0\n"
complete_tests = "When this option is choosen the program will undergo\n" \
                 " all the tests implemented in the tests modules."
mutants_from_file_description = "This option specifies a path to read the output names and the\n" \
                                " mutations from a file.\n" \
                                "The file should have one mutation per line and no blank lines.\n" \
                                "Each line should have the following format:\n" \
                                "  output_name-->XXX N C YYY\n" \
                                "    --> indicates a tab space\n" \
                                "    XXX is the initial aminoacid in 3 letters code\n" \
                                "    N the residue number\n" \
                                "    C the chain to modify\n" \
                                "    YYY the desired aminoacid in 3 letters code\n"
mutant_multiple_description = "This parameter creates one structure with all the\n" \
                              " mutations specified. But it doesn't check for clashes."
charge_terminals_description = "When this option is present the program will charge all the terminals residues," \
                               "including those in gaps."
remove_terminal_missing_description = "Whe this option is chosen the program will remove the terminal residues " \
                                      "in case they are missing heavy atoms in the backbone of the protein" \
                                      "(that is: N, CA, C, O)."
no_gaps_ter_description = "When this option is present the program won't add a TER mark whenever it finds " \
                       "a gap in the sequence."
pdb_resolution_description = "A float indicating the pdb resolution, this number will be used as the maximum allowed " \
                             "between the N atom of a residue and the C atom of the previous residue when checking" \
                             "for gaps. In this check the minimum distance that will be used is 1.55A."
make_unique_description = "This option is optional and should be present if the user wants to make the names " \
                          "in a given CHAIN unique."
remove_terminal_missing = "When this option is present the program will eliminate the first or last residue if they " \
                          "have missing atoms, and fix the next or previous residues to be the terminal one."