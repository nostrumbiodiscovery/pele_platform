#!/bin/python
import sys
from prody import *

from pele_platform.Utilities.PPP.adjustments_module import WritingAtomNames, FixStructureResnames, FixAtomNames, SolveClashes
from pele_platform.Utilities.PPP.checks_module import CheckMutation, CheckClashes
from pele_platform.Utilities.PPP.checks_module import CheckStructure, CheckforGaps
from pele_platform.Utilities.PPP.global_processes import ParseArguments, FindInitialAndFinalResidues, PDBwriter, RenumberStructure
from pele_platform.Utilities.PPP.hydrogens_addition import FixStructure
from pele_platform.Utilities.PPP.mutational_module import Mutate
from pele_platform.Utilities.PPP.global_variables import coordination_geometries

__author__ = 'jelisa'

"""
This block adds the hid, hip, and hie residues to prody. Otherwise this would consider
this aminoacids as heteroatoms.
"""
addNonstdAminoacid('HID', 'aromatic', 'basic', 'cyclic', 'large', 'polar', 'surface')
addNonstdAminoacid('HIE', 'aromatic', 'basic', 'cyclic', 'large', 'polar', 'surface')
addNonstdAminoacid('HIP', 'aromatic', 'basic', 'cyclic', 'large', 'polar', 'surface')
addNonstdAminoacid('CYT', 'neutral', 'acyclic', 'medium', 'polar', 'buried')
addNonstdAminoacid('LYN', 'neutral', 'acyclic', 'large', 'polar', 'buried')


def main(input_pdb, pdb_resolution, output_pdb="", no_gaps_ter=False, charge_terminals=False, make_unique=False,
         remove_terminal_missing=False, mutant_multiple=False, mutation=""):
    try:
        initial_structure = parsePDB(input_pdb)
    except IOError:
        print("The file '{}' isn't a valid file\nCheck that it does exist and try again.".format(input_pdb))
        sys.exit()
    initial_residue, final_residue = FindInitialAndFinalResidues(initial_structure)
    # ff_parameters = ReadForceFieldParameters(force_field)

    print("* Checking for insertion codes.")
    insertion_codes = [icode for icode in initial_structure.getIcodes() if icode]
    if insertion_codes:
        print("  * Due to the insertion codes the structure will be RENUMBERED starting from 1 for each chain.\n" \
              "  ** Otherwise PELE would fail.")
        structure2use = RenumberStructure(initial_structure)
    else:
        structure2use = initial_structure
    print("* Checking for gaps.")
    gaps, not_gaps = CheckforGaps(structure2use, pdb_resolution)
    if gaps is None and not_gaps is None:
        print("WARNING: Problems when checking for gaps, so don't trust the existence of gaps.")
        gaps, not_gaps = {}, {}
    print("* Checking and Fixing the Residues Names:")
    structure2use = FixStructureResnames(structure2use, make_unique)
    print("* Checking and fixing the Atoms Names:")
    structure2use = FixAtomNames(structure2use, gaps, not_gaps)
    print("* Checking the structure for missing atoms:")
    residues2fix, residues2remove, metals2coordinate, residues_without_template = CheckStructure(structure2use, gaps,
                                                                                                 not_gaps,
                                                                                                 charge_terminals,
                                                                                                 remove_terminal_missing)
    if residues2fix:
        print('* Placing the missing atoms and removing the extra atoms:')
        structure2use = FixStructure(structure2use, residues2fix, gaps, charge_terminals)

    if not mutation:
        print('Writing the structure to {}'.format(output_pdb[0]))
        if make_unique:
            ligand_chain = structure2use.select("chain {}".format(make_unique))
            if ligand_chain:
                not_proteic_ligand = structure2use.select("chain {}".format(make_unique)).hetero
            else:
                not_proteic_ligand = None
            PDBwriter(output_pdb[0], WritingAtomNames(structure2use), make_unique, residues2remove,
                      no_gaps_ter, not_proteic_ligand, gaps, not_gaps)
        else:
            not_proteic_ligand = None
            PDBwriter(output_pdb[0], WritingAtomNames(structure2use), make_unique, residues2remove,
                      no_gaps_ter, not_proteic_ligand, gaps, not_gaps)

        return residues_without_template, gaps, metals2coordinate

if __name__ == '__main__':
    arguments = ParseArguments()
    if arguments is None:
        sys.exit()
    else:
        main(arguments.input_pdb, arguments.pdb_resolution, output_pdb=arguments.output_pdb,
             no_gaps_ter=arguments.no_gaps_ter, charge_terminals=arguments.charge_terminals,
             make_unique=arguments.make_unique, remove_terminal_missing=arguments.remove_terminal_missing,
mutant_multiple=arguments.mutant_multiple, mutation=arguments.mutation)
