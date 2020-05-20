import pele_platform.constants.constants as cs
from Bio.PDB import PDBParser, NeighborSearch, Selection, Vector, vectors, PDBIO
import itertools
import numpy as np
import re

def find_metal_geo(protein_file):

    # read in the protein file
    parser = PDBParser()
    structure = parser.get_structure("protein", protein_file)
    structure_list =  Selection.unfold_entities(structure, "A")

    # find metals
    metals = []
    for chain in structure.get_chains():
        for residue in chain.get_residues():
            for atom in residue.get_atoms():
                if atom.element in cs.metals:
                    metals.append([atom, residue, chain])

    # check metal contacts
    coordinated_atoms = []
    chains = []
    residues = []
    angles = []
    atoms = []
    contacts = []
    contacts_ang = []

    for metal in metals:
        coords = metal[0].coord
        contacts_chain = NeighborSearch(structure_list).search(coords, 4.0, "C")
        
        for chain in contacts_chain:
            contacts_residue = NeighborSearch(structure_list).search(coords, 4.0, "R")
            
            for residue in contacts_residue:
                contacts_atoms = NeighborSearch(structure_list).search(coords, 4.0, "A")
                contacts_atoms = [c for c in contacts_atoms if c.element not in [metal[0].name, "C", "H"]] # exclude self-contacts, carbons and hydrogens
                
                for atom in contacts_atoms:
                    if residue in chain.get_residues() and atom in residue.get_atoms():
                        contacts.append([atom, residue, chain])
       
        combinations = list(itertools.combinations(contacts, 2))
        combinations = [list(c) for c in combinations] 
        
        # get all atom - metal - atom angles
        for c in combinations:
            vi = Vector(c[0][0].coord)
            vj = Vector(c[1][0].coord)
            angle = vectors.calc_angle(vi, coords, vj) * 180 / np.pi
            c.append(angle)
    
        # angle classification
        ang_180 = []
        ang_90 = []
        ang_109 = []
        lower = 0.9
        upper = 1.1

        for c in combinations:
            a = c[2]
            if 180*lower<=a<=180*upper:
                ang_180.append(c)
            if 90*lower<=a<=90*upper:
                ang_90.append(c)
            if 109.5*lower<=a<=109.5*upper:
                ang_109.append(c)
    
        # check geometries
        if len(ang_180) >= 3 and len(ang_90) >=12:
            geo = "octahedral"
            coordinated_atoms.extend(ang_180)
            coordinated_atoms.extend(ang_90)
        elif len(ang_180) >= 2 and len(ang_90) >=4:
            geo = "square planar"
            coordinated_atoms.extend(ang_180)
            coordinated_atoms.extend(ang_90)
        elif len(ang_109) >= 6:
            geo = "tetrahedral"
            coordinated_atoms.extend(ang_190)
        else:
            raise Warning("Failed to determine geometry around the metal atom. Add constraints manually.")
        
        print("Found {} geometry around the metal center. Adding constraints.".format(geo))
        
        # format string
        yaml_string = "{}-{}-{}:{}:{}-{}:{}:{}"
        output = []
        spring_const = 50

        for c in coordinated_atoms:
            atom1, atom2, angle = c
            atomname1 = atom1[0].name
            resnum1 = atom1[1].get_id()[1]
            chain1 = atom1[2].get_id()

            atomname2 = metal[0].name
            resnum2 = metal[1].get_id()[1]
            chain2 = metal[2].get_id()

            atom_dist = atom1[0] - metal[0]
            out = yaml_string.format(spring_const, atom_dist, chain1, resnum1, atomname1, chain2, resnum2, atomname2)
            output.append(out)
        output = list(set(output))
        
        for o in output:
            print(o)
        
    return output


def main(protein_file):
    find_metal_geo(protein_file)

if __name__ == "__main__":
    main("1zop_angles.pdb")
