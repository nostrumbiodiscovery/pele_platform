import sys
from prody import *
import os
from adjustments_module import WritingAtomNames, FixStructureResnames, FixAtomNames, SolveClashes
from checks_module import CheckMutation, CheckClashes
from checks_module import CheckStructure, CheckforGaps
from global_processes import ParseArguments, FindInitialAndFinalResidues, PDBwriter, RenumberStructure
from hydrogens_addition import FixStructure
from mutational_module import Mutate

__author__ = 'jelisa'




def main(input_pdb, output_pdb = [], mutation = False, no_gaps_ter=False, make_unique=False, remove_terminal_missing=False):


    """
    This block adds the hid, hip, and hie residues to prody. Otherwise this would consider
    this aminoacids as heteroatoms.
    """
    addNonstdAminoacid('HID', 'aromatic', 'basic', 'cyclic', 'large', 'polar', 'surface')
    addNonstdAminoacid('HIE', 'aromatic', 'basic', 'cyclic', 'large', 'polar', 'surface')
    addNonstdAminoacid('HIP', 'aromatic', 'basic', 'cyclic', 'large', 'polar', 'surface')
    addNonstdAminoacid('CYT', 'neutral', 'acyclic', 'medium', 'polar', 'buried')
    addNonstdAminoacid('LYN', 'neutral', 'acyclic', 'large', 'polar', 'buried')


    if not output_pdb:
        output_pdb.append("{}_processed.pdb".format(os.path.splitext(input_pdb)[0]))
  

    if input_pdb is None:
        sys.exit()

    try:
        initial_structure = parsePDB(input_pdb)
    except IOError:
        print "The file '{}' isn't a valid file\nCheck that it does exist and try again.".format(input_pdb)
        sys.exit()
    initial_residue, final_residue = FindInitialAndFinalResidues(initial_structure)
    # ff_parameters = ReadForceFieldParameters(force_field)

    print "* Checking for gaps."
    gaps, not_gaps = CheckforGaps(initial_structure)
    if gaps is None and not_gaps is None:
        print "WARNING: Problems when checking for gaps, so don't trust the existence of gaps."
        gaps, not_gaps = {}, {}
    print "* Checking for insertion codes."
    insertion_codes = [icode for icode in initial_structure.getIcodes() if icode]
    if insertion_codes:
        print " *The structure will be renumbered starting from 1 for each chain."
        structure2use = RenumberStructure(initial_structure, gaps, not_gaps)
    else:
        structure2use = initial_structure
    print "* Checking and Fixing the Residues Names:"
    structure2use = FixStructureResnames(structure2use, make_unique)
    print "* Checking and fixing the Atoms Names:"
    structure2use = FixAtomNames(structure2use, gaps, not_gaps)
    print "* Checking the structure for missing atoms:"
    residues2fix, residues2remove, missing_res = CheckStructure(structure2use, gaps, not_gaps, remove_terminal_missing)
    if residues2fix:
        print '* Placing the missing atoms:'
        structure2use = FixStructure(structure2use, residues2fix)
    print mutation

    if not mutation:
        print 'Writing the structure to {}'.format(output_pdb[0])
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

    return os.path.abspath(output_pdb[0]), missing_res 

if __name__ == '__main__':
   # arguments = ParseArguments()
    main(sys.argv[1])
