#!/bin/python
import os
import sys
from prody import *

from pele_platform.Utilities.PPP.adjustments_module import WritingAtomNames, FixStructureResnames, FixAtomNames, SolveClashes
from pele_platform.Utilities.PPP.checks_module import CheckMutation, CheckClashes
from pele_platform.Utilities.PPP.checks_module import CheckStructure, CheckforGaps
from pele_platform.Utilities.PPP.global_processes import ParseArguments, FindInitialAndFinalResidues, PDBwriter, RenumberStructure
from pele_platform.Utilities.PPP.hydrogens_addition import FixStructure
from pele_platform.Utilities.PPP.mutational_module import Mutate
from pele_platform.Utilities.PPP.global_variables import coordination_geometries
import pele_platform.Utilities.Helpers.constraints as ct

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


def main(input_pdb, pele_dir, output_pdb=["",],  no_gaps_ter=False, charge_terminals=False, make_unique=False,
         remove_terminal_missing=False, mutant_multiple=False, mutation="", mid_chain_nonstd_residue=[], skip=False, resolution=2.5):

    if not output_pdb[0]:
        output = os.path.splitext(os.path.basename(input_pdb))[0]
        output_pdb[0] = os.path.join(pele_dir,"{}_processed.pdb".format(output))

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
    gaps, not_gaps = CheckforGaps(structure2use, resolution)
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
            if not skip:
                PDBwriter(output_pdb[0], WritingAtomNames(structure2use), make_unique, residues2remove,
                      no_gaps_ter, not_proteic_ligand, gaps, not_gaps, mid_chain_nonstd_residue)
            else:
                shutil.copy(input_pdb, output_pdb[0])

        ########METAL ATOMS PROCESS###############
        coordinated_atoms_ids = {}
        for metal, atoms_list in metals2coordinate.items():
            metal_id = "{} {} {}".format(WritingAtomNames(metal).getNames()[0].replace(' ','_'),
                                         metal.getChid(),
                                         metal.getResnum())
            atoms_ids = [["{} {} {} {}".format(at.getResname(),
                                               at.getResnum(), at.getChid(),
                                               WritingAtomNames(at).getNames()[0].replace(' ', '_'),),
                          calcDistance(metal, at)[0]] for at in atoms_list]
            if len(atoms_list) in [x[1] for x in coordination_geometries.itervalues()]:
                coordinated_atoms_ids[metal_id] = atoms_ids

        ###########RETRIEVE CONSTANTS###############
        constr = ct.retrieve_constraints(output_pdb[0], gaps, coordinated_atoms_ids)



        return output_pdb[0], residues_without_template, gaps, coordinated_atoms_ids, constr

if __name__ == '__main__':
    arguments = ParseArguments()
    if arguments is None:
        sys.exit()
    else:
        main(arguments.input_pdb, arguments.pdb_resolution, output_pdb=arguments.output_pdb,
             no_gaps_ter=arguments.no_gaps_ter, charge_terminals=arguments.charge_terminals,
             make_unique=arguments.make_unique, remove_terminal_missing=arguments.remove_terminal_missing,
mutant_multiple=arguments.mutant_multiple, mutation=arguments.mutation)
