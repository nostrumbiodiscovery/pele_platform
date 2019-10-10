#!/bin/python
import sys
from prody import *
import os

import msm_pele.Helpers.constraints as ct
from msm_pele.PPP.adjustments_module import WritingAtomNames, FixStructureResnames, FixAtomNames, SolveClashes
from msm_pele.PPP.checks_module import CheckMutation, CheckClashes
from msm_pele.PPP.checks_module import CheckStructure, CheckforGaps
from msm_pele.PPP.global_processes import ParseArguments, FindInitialAndFinalResidues, PDBwriter, RenumberStructure
from msm_pele.PPP.hydrogens_addition import FixStructure
from msm_pele.PPP.mutational_module import Mutate
from msm_pele.PPP.global_variables import coordination_geometries

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


def main(input_pdb, pele_dir, output_pdb=["",], no_gaps_ter=False, charge_terminals=False, make_unique=False,
         remove_terminal_missing=False, mutant_multiple=False, mutation="", mid_chain_nonstd_residue=[], renumber=True, dynamic_waters=[]):
    if not output_pdb[0]:
        output = os.path.splitext(os.path.basename(input_pdb))[0]
        output_pdb[0] = os.path.join(pele_dir,"{}_processed.pdb".format(output))
        print(output_pdb[0])
    try:
        initial_structure = parsePDB(input_pdb)
    except IOError:
        print("The file '{}' isn't a valid file\nCheck that it does exist and try again.".format(input_pdb))
        sys.exit()
    initial_residue, final_residue = FindInitialAndFinalResidues(initial_structure)
    # ff_parameters = ReadForceFieldParameters(force_field)

    print("* Checking for gaps.")
    gaps, not_gaps = CheckforGaps(initial_structure)
    if gaps is None and not_gaps is None:
        print("WARNING: Problems when checking for gaps, so don't trust the existence of gaps.")
        gaps, not_gaps = {}, {}
    print("* Checking for insertion codes.")
    insertion_codes = [icode for icode in initial_structure.getIcodes() if icode]
    if insertion_codes and renumber:
        print(" *The structure will be renumbered starting from 1 for each chain.")
        structure2use = RenumberStructure(initial_structure, gaps, not_gaps)
    else:
        structure2use = initial_structure
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
    print(mutation)

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
                      no_gaps_ter, not_proteic_ligand, gaps, not_gaps, mid_chain_nonstd_residue)

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
        constr = ct.retrieve_constraints(output_pdb[0], gaps, coordinated_atoms_ids, dynamic_waters=dynamic_waters)
        return output_pdb[0], residues_without_template, gaps, coordinated_atoms_ids, constr
    else:
        clashes = []
        mutated_structure = None
        for mutation, output_file in zip(mutation, output_pdb):
            print('* Checking the mutation:')
            print(" Mutation: {0[ini_resname]} {0[chain]} {0[resnum]} {0[fin_resname]}".format(mutation))
            correct_mutation = CheckMutation(structure2use, mutation)
            if not correct_mutation:
                exit_message = "The mutation was incorrect, check your parameters.\n" \
                               "The checked structure will be written to {}".format(output_file)
                PDBwriter(output_file, WritingAtomNames(structure2use), make_unique, gaps, no_gaps_ter, not_gaps)
                continue
            else:
                print("Output_file name: {0}".format(output_file))
                mutated_structure, zmatrix = Mutate(structure2use, mutation)
                if not mutant_multiple:
                    if mutation[0]['fin_resname'] in ["ALA", "GLY"]:
                        print("The ALA and the GLY don't have any rotamer to try.")
                    else:
                        print("Checking Clashes:")
                        try:
                            clashes = CheckClashes(mutated_structure, mutation, zmatrix,
                                                   initial_residue, final_residue)
                        except ValueError:
                            pass
                        else:
                            if not clashes:
                                print("Structure without clashes.")
                            else:
                                mutated_structure = SolveClashes(mutated_structure, clashes,
                                                                 mutation, zmatrix,
                                                                 initial_residue, final_residue)
                    mutated_structure.setTitle("mutated structure")
                    PDBwriter(output_file, WritingAtomNames(mutated_structure), make_unique, gaps, no_gaps_ter,
                              not_gaps)
                else:
                    print("Multiple mutations at the same time are still under development.")
                    structure2use = mutated_structure
        if mutant_multiple and mutated_structure is not None:
            PDBwriter(output_pdb, WritingAtomNames(mutated_structure), gaps, not_gaps, no_gaps_ter)
    # return residues_without_template, gaps, metals2coordinate


if __name__ == '__main__':
    arguments = ParseArguments()
    if arguments is None:
        sys.exit()
    else:
        main(arguments.input_pdb, output_pdb=arguments.output_pdb, no_gaps_ter=arguments.no_gaps_ter,
             charge_terminals=arguments.charge_terminals, make_unique=arguments.make_unique,
             remove_terminal_missing=arguments.remove_terminal_missing, mutant_multiple=arguments.mutant_multiple,
             mutation=arguments.mutation)
