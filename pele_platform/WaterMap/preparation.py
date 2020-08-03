import numpy as np
import os
from Bio.PDB import PDBParser, PDBIO, NeighborSearch, Selection
from pele_platform.Utilities.Helpers.randomize import calculate_com


def prepare_system(protein_file, user_center, user_radius):

    pdb_no_waters = remove_water(protein_file)
    pdb_with_water = add_water_box(pdb_no_waters, user_center)
    equilibration_input = remove_overlaps(pdb_with_water, user_center, user_radius)

    return equilibration_input


def remove_water(protein_file):

    new_lines = []
    new_protein_file = os.path.basename(protein_file).replace(".pdb", "_prep.pdb")
    new_protein_file = os.path.abspath(new_protein_file)

    with open(protein_file, "r") as file:
        lines = file.readlines()

        for line in lines:
            if line[17:20] != "HOH":
                new_lines.append(line)

    with open(new_protein_file, "w") as new_file:

        for line in new_lines:
            new_file.write(line)

    return new_file.name


def add_water_box(protein_file, user_center):

    water_box_file = "water_box.pdb"  # make more general
    output_file = "translated_water.pdb"

    parser = PDBParser()
    water_box = parser.get_structure("water_box", water_box_file)

    # calculate translation vector
    water_box_center = calculate_com(water_box)
    move_vector = water_box_center - user_center

    # translate water box
    for atom in water_box.get_atoms():
        new_position = np.array(list(atom.get_vector())) - move_vector
        atom.set_coord(new_position)

    # save new water box coordinates
    io = PDBIO()
    io.set_structure(water_box)
    temp_output_file = "translated_water_temp.pdb"
    io.save(temp_output_file)

    # merge translated water file with the protein PDB
    to_save = []

    with open(protein_file, "r") as protein_file:
        lines = protein_file.readlines()
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("TER"):
                to_save.append(line)

    with open(temp_output_file, "r") as water_file:
        lines = water_file.readlines()
        for line in lines:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                to_save.append(line)

    with open(output_file, "w+") as out_file:
        for i in to_save:
            out_file.write(i)
        out_file.write("\nTER\nEND")

    return output_file


def remove_overlaps(protein_file, user_center, user_radius):

    equilibration_input = "equilibration_input.pdb"

    parser = PDBParser()
    structure = parser.get_structure("structure", protein_file)
    protein_list = Selection.unfold_entities(structure, "A")
    residues_to_remove = []
    to_remove = []

    # check if water residues are inside user-defined sphere and are not overlapping with protein
    for residue in structure.get_residues():
        if residue.resname == "TIP":
            for atom in residue.get_atoms():

                if atom.name == "OH2":
                    oxygen_coords = list(atom.get_vector())
                    if np.linalg.norm(np.subtract(oxygen_coords, user_center)) > user_radius:
                        residues_to_remove.append(residue.get_id()[1])
                    neighbours = NeighborSearch(protein_list).search(oxygen_coords, 1.0, "R")

                    for contact in neighbours:
                        if contact.get_resname() != "TIP":
                            residues_to_remove.append(residue.get_id()[1])

    residues_to_remove = list(set(residues_to_remove))  # removing duplicates
    residues_to_remove = [str(a) for a in residues_to_remove]

    with open(protein_file, "r") as file:
        lines = file.readlines()
        for line in lines:
            if line[22:26].strip() in residues_to_remove and len(line[20:22].strip()) == 0:
                to_remove.append(line)

    final_lines = [line for line in lines if line not in to_remove]

    with open(equilibration_input, "w+") as file_out:
        for line in final_lines:
            file_out.write(line)

    return equilibration_input


if __name__ == "__main__":
    prepare_system("fake_input.pdb", [29.8419990540, 48.9189987183, 61.5859985352], 6.0)
