import pele_platform.constants.constants as cs
from Bio.PDB import PDBParser, NeighborSearch, Selection, Vector, vectors, PDBIO
import itertools
import numpy as np
import re
import warnings

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
    output = []

    for metal in metals:
        coords = metal[0].coord
        coordinated_atoms = []
        contacts = []

        for chain in structure.get_chains():
            
            for residue in chain.get_residues():
                contacts_atoms = NeighborSearch(structure_list).search(coords, 3.2, "A")
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
        lower = 0.80
        upper = 1.20

        for c in combinations:
            a = c[2]
            if 180*lower<=a<=180*upper:
                ang_180.append(c)
            if 90*lower<=a<=90*upper:
                ang_90.append(c)
            if 109.5*lower<=a<=109.5*upper:
                ang_109.append(c)
        
        # check geometries
        if len(ang_180) == 3 and len(ang_90) == 12:
            geo = "octahedral"
            coordinated_atoms.extend(ang_180)
            coordinated_atoms.extend(ang_90)
        elif len(ang_180) == 2 and len(ang_90) == 4:
            geo = "square planar"
            coordinated_atoms.extend(ang_180)
            coordinated_atoms.extend(ang_90)
        elif len(ang_109) == 6:
            geo = "tetrahedral"
            coordinated_atoms.extend(ang_109)
        else:
            geo = None
            warnings.warn(("Failed to determine geometry around {} (residue {}). Add constraints manually.".format(metal[0].name, metal[1].get_id()[1])), Warning)
       
        if geo:
            print("Found {} geometry around {} (residue {}). Adding constraints.".format(geo, metal[0].name, metal[1].get_id()[1]))
        
        # format string
        yaml_string = "{}-{}-{}:{}:{}-{}:{}:{}"
        spring_const = 50

        string_atoms = []
        for c in coordinated_atoms:
            atom1, atom2, angle = c

            if atom1 not in string_atoms:
                string_atoms.append(atom1)
            if atom2 not in string_atoms:
                string_atoms.append(atom2)

        for atom in string_atoms:
            atomname1 = atom[0].name
            resnum1 = atom[1].get_id()[1]
            chain1 = atom[2].get_id()

            atomname2 = metal[0].name
            resnum2 = metal[1].get_id()[1]
            chain2 = metal[2].get_id()

            atom_dist = atom[0] - metal[0]
            out = yaml_string.format(spring_const, atom_dist, chain1, resnum1, atomname1, chain2, resnum2, atomname2)
            output.append(out)
        output = list(set(output))
        if output:
            output = ['{}'.format(o) for o in output]
            #output = ", ".join(output)
    return output


def main(protein_file):
    output = find_metal_geo(protein_file)
