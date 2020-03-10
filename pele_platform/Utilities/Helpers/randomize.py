from Bio.PDB import PDBParser, PDBIO, NeighborSearch, Selection
from Bio.PDB.vectors import rotaxis2m, Vector
import numpy as np
import os


def randomize_starting_position(ligand_file, complex_file, nposes=200):
    """
    Randomize initial ligand position around the receptor.
    Default number of poses = 200.
    :param ligand_file:
    :param complex_file:
    :param nposes:
    :return:
    """

    parser = PDBParser()
    output = []
    structure = parser.get_structure('protein', complex_file)
    ligand = parser.get_structure('ligand', ligand_file)

    # calculate protein COM
    structure_mass = 0.0
    com_protein = np.zeros(3)
    for atom in structure.get_atoms():
        com_protein = com_protein + np.array(list(atom.get_vector())) * atom.mass
        structure_mass += atom.mass
    com_protein = com_protein / structure_mass

    # calculate ligand COM
    ligand_mass = 0.0
    com_ligand = np.zeros(3)
    for atom in ligand.get_atoms():
        com_ligand = com_ligand + np.array(list(atom.get_vector())) * atom.mass
        ligand_mass += atom.mass
    com_ligand = com_ligand / ligand_mass

    # calculate vector to move the ligand
    move_vector = com_ligand - com_protein

    # translate the ligand to the protein COM
    original_coords = []
    for atom in ligand.get_atoms():
        ligand_origin = np.array(list(atom.get_vector())) - move_vector
        original_coords.append(ligand_origin)
        atom.set_coord(ligand_origin)

    # calculating the maximum radius of the protein from the origin
    coor = []
    for atom in structure.get_atoms():
        coor.append(list(atom.get_vector() - com_protein))

    coor = np.array(coor)
    coor_max = np.amax(coor, axis=0)

    d = np.sqrt(np.sum(coor_max ** 2))

    D = np.ceil(6.0 + d)  # radius of the sphere from the origin
    print("Sampling {}A spherical box around the centre of the receptor".format(D))
    sphere_cent = com_protein

    j = 0
    print("Generating {} poses...".format(nposes))
    while (j < nposes):
        phi = np.random.uniform(0, 2 * np.pi)
        costheta = np.random.uniform(-1, 1)
        u = np.random.uniform(0, 1)
        theta = np.arccos(costheta)

        r = D * np.cbrt(u)
        x = r * np.sin(theta) * np.cos(phi)
        y = r * np.sin(theta) * np.sin(phi)
        z = r * np.cos(theta)

        # move ligand to the starting point (protein COM)
        for atom, coord in zip(ligand.get_atoms(), original_coords):
            atom.set_coord(coord)

        translation = (x, y, z)

        for atom in ligand.get_atoms():
            new_pos_lig_trans = np.array(list(atom.get_vector())) - translation
            atom.set_coord(new_pos_lig_trans)

        new_ligand_COM = np.zeros(3)
        ligand_mass = 0
        for atom in ligand.get_atoms():
            new_ligand_COM = new_ligand_COM + np.array(list(atom.get_vector())) * atom.mass
            ligand_mass += atom.mass
        new_ligand_COM = new_ligand_COM / ligand_mass

        # check if it's inside the sampling sphere
        dist = np.sqrt((new_ligand_COM[0] - sphere_cent[0]) ** 2 + (new_ligand_COM[1] - sphere_cent[1]) ** 2 + (
                new_ligand_COM[2] - sphere_cent[2]) ** 2)

        if dist < D:
            angles = (np.random.uniform(0, 360), np.random.uniform(0, 360), np.random.uniform(0, 360))

            # random rotation
            translation = np.zeros(3)  # blank
            for atom in ligand.get_atoms():
                rotation1 = rotaxis2m(angles[0], Vector(1, 0, 0))
                atom.transform(rotation1, translation)
                rotation2 = rotaxis2m(angles[1], Vector(0, 1, 0))
                atom.transform(rotation2, translation)
                rotation3 = rotaxis2m(angles[2], Vector(0, 0, 1))
                atom.transform(rotation3, translation)
                final_coords = np.array(list(atom.get_vector()))

            # check contacts at: 5A (no contacts) and 8A (needs contacts)
            contacts5 = []
            protein_list = Selection.unfold_entities(structure, "A")

            for atom in ligand.get_atoms():
                center = np.array(list(atom.get_vector()))
                contacts5 = NeighborSearch(protein_list).search(center, 5.0, "A")
                contacts5_results = []
                if contacts5:
                    contacts5_results.append(False)
                    print("Contacts at 5A.")
                else:
                    contacts5_results.append(True)

            contacts8 = []
            contacts8_results = []

            for atom in ligand.get_atoms():
                center = np.array(list(atom.get_vector()))
                contacts8 = NeighborSearch(protein_list).search(center, 8.0, "A")
                if contacts8:
                    contacts8_results.append(True)
                else:
                    contacts8_results.append(False)
                    print("No contacts at 8A.")
            if (any(contacts8_results) and all(contacts5_results)):
                j += 1
                io = PDBIO()
                io.set_structure(ligand)
                output.append(io.save('correct_ligand{}.pdb'.format(j)))
        else:
            print("Ligand it not inside the sampling sphere.")
    print("{} poses created successfully.".format(j))
    return output
