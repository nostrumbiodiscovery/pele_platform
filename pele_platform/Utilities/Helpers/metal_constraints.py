import pele_platform.constants.constants as cs
from Bio.PDB import PDBParser, NeighborSearch, Selection, Vector, vectors, PDBIO
import itertools
import numpy as np
import re

def find_metal_geo(protein_file, permissive=False, external=None):

    checked_metals = []

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

        metal_str = "{}:{}:{}".format(metal[2], metal[1].get_id(), metal[0])

        if metal_str not in external and list(metal[0].coord) not in checked_metals:
            checked_metals.append(list(metal[0].coord))
            print("checked_metals", checked_metals)
            coords = metal[0].coord
            coordinated_atoms = []
            contacts = []

            for chain in structure.get_chains():
                
                for residue in chain.get_residues():
                    contacts_atoms = NeighborSearch(structure_list).search(coords, 2.6, "A")
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
            
            if permissive:
                lower = 0.65
                upper = 1.35
            else:
                lower = 0.8
                upper = 1.2

            for c in combinations:
                a = c[2]
                print("c",c)
                if 180*lower<=a<=180*upper:
                    ang_180.append(c)
                    print("180")
                if 90*lower<=a<=90*upper:
                    ang_90.append(c)
                    print("90")
                if 109.5*lower<=a<=109.5*upper:
                    ang_109.append(c)
                    print("109.5")
            print("180", len(ang_180), "109.5", len(ang_109), "90", len(ang_90))

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
                if not permissive:
                    raise Exception("Failed to determine geometry around {} (residue {}). Set 'permissive_metal_constr: true' to allow more permissive angle classification or add constraints manually.".format(metal[0].name, metal[1].get_id()[1]))
                else:
                    raise Exception("Failed to determine geometry around {} (residue {}). Add constraints manually.".format(metal[0].name, metal[1].get_id()[1]))
           
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
    
    return output
