from Bio.PDB import PDBParser, PDBIO, NeighborSearch, Selection
from Bio.PDB.vectors import Vector
import numpy as np
import os, argparse

def randomize_starting_position(ligand_file, complex_file, outputfolder=".", nposes=200, test=False):
    """
    Randomize initial ligand position around the receptor.
    Default number of poses = 200.
    :param ligand_file:
    :param complex_file:
    :param nposes:
    :return:
    """
    if test:  np.random.seed(42)

    # read in files

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

    # calculating the maximum d of the ligand
    coor_ligand = []
    for atom in ligand.get_atoms():
        coor_ligand.append(list(atom.get_vector() - com_ligand))

    coor_ligand = np.array(coor_ligand)
    coor_ligand_max = np.amax(coor_ligand, axis=0)

    d_ligand = np.sqrt(np.sum(coor_ligand_max ** 2))

    # set threshold for near and far contacts based on ligand d
    if d_ligand < 5.0:
        d5_ligand = 5.0
    else:
        d5_ligand = d_ligand / 2

    if d_ligand > 8.0:
        d8_ligand = d_ligand / 2 + 3
    else:
        d8_ligand = 8.0

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

    # radius of the sphere from the origin
    D = np.ceil(6.0 + d)
    print("Sampling {}A spherical box around the centre of the receptor".format(D))
    sphere_cent = com_protein

    j = 0
    print("Generating {} poses...".format(nposes))
    while (j < nposes):

        # generate random coordinates
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

        # translate ligand to a random position
        translation = (x, y, z)
        for atom in ligand.get_atoms():
            new_pos_lig_trans = np.array(list(atom.get_vector())) - translation
            atom.set_coord(new_pos_lig_trans)

        # calculate ligand COM in the new position
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

            # check contacts at: 5A (no contacts) and 8A (needs contacts)
            protein_list = Selection.unfold_entities(structure, "A")
            contacts5 = []
            contacts5 = NeighborSearch(protein_list).search(new_ligand_COM, d5_ligand, "A")
            contacts8 = []
            contacts8 = NeighborSearch(protein_list).search(new_ligand_COM, d8_ligand, "A")

            if contacts8 and not contacts5:
                j += 1
                io = PDBIO()
                io.set_structure(ligand)
                output_name = os.path.join(outputfolder, 'ligand{}.pdb'.format(j))
                io.save(output_name)
                output.append(output_name)
    print("{} poses created successfully.".format(j))
    return output,  D, list(sphere_cent)


def join(receptor, ligands, residue, output_folder=".", output="input{}.pdb"):
    """
    Join receptor and ligand PDB files into one, conserving old formatting
    and not repeating atom numbers.
    """

    with open(receptor, "r") as f:
        lines = f.readlines()
        receptor_content = [line for line in lines if line[17:20] != residue and not line.startswith("CONECT") and not line.startswith("END")]
        connects = [line for line in lines if line.startswith("CONECT")]
        ligand_content_without_coords = [line[0:27] + "{}" + line[56:] for line in lines if line[17:20] == residue]
        atom_nums = [line[6:11] for line in lines if line[17:20] == residue]

    outputs = []
    for i, ligand in enumerate(ligands):
        with open(ligand, "r") as fin:
            # exclude connects but keep initial atom names (CL problem)
            ligand_coords = {line[12:16].strip(): line[27:56] for line in fin if
                             line.startswith("ATOM") or line.startswith("HETATM")}
            #assert len(ligand_coords) == len(ligand_content_without_coords), "Experimental part - send an issue to github"

            ligand_content = []
            for pdb_block in ligand_content_without_coords:
                atom_name = pdb_block[12:16].strip()
                coord = ligand_coords[atom_name]
                ligand_pdb_line = pdb_block.format(coord)
                ligand_content.append(ligand_pdb_line)


            for j, line in enumerate(ligand_content):
                ligand_content[j] = line[:6] + "{:>5}".format(atom_nums[j]) + line[11:]
        content_join_file = receptor_content + ["TER\n"] + ligand_content + ["TER\n"] + connects + ["END"]
        output_path = os.path.join(output_folder, output.format(i))
        with open(output_path, "w") as fout:
            fout.write("".join(content_join_file))
        outputs.append(output_path)

    return outputs


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--ligand", type=str, required=True, help="Ligand pdb file")
    parser.add_argument("--receptor", type=str, required=True, help="Receptor pdb file")
    parser.add_argument("--resname", type=str, required=True, help="Ligand resname")
    parser.add_argument("--poses", type=int, default=20, help="How many input poses to produce")
    parser.add_argument("--output_folder", type=str, default=".", help="output folder")
    args = parser.parse_args()
    return os.path.abspath(args.ligand), os.path.abspath(args.receptor), args.resname, args.poses, args.output_folder


if __name__ == "__main__":
    ligand, receptor, resname, poses, output_folder = parse_args()
    output, sphere_cent, D = randomize_starting_position(ligand, "input_ligand.pdb", resname, receptor, None, None,
                                                         output_folder, poses=poses)
    # Make format ready to be captured"
    print("OUTPUT; {}; {}; {}".format(output, sphere_cent, D))
