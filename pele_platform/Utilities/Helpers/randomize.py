from Bio.PDB import PDBParser, PDBIO, NeighborSearch, Selection, rotaxis
from Bio.PDB.vectors import Vector
import numpy as np
import os, argparse
import time
import pele_platform.Errors.custom_errors as cs

def calculate_com(structure):
    structure_mass = 0.0
    com = np.zeros(3)
    for atom in structure.get_atoms():
        com = com + np.array(list(atom.get_vector())) * atom.mass
        structure_mass += atom.mass
    com = com / structure_mass
    return com


def randomize_starting_position(ligand_file, complex_file, outputfolder=".", nposes=200, test=False, user_center=None,
                                logger=None):
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
    COI = np.zeros(3)

    # get center of interface (if PPI)
    if user_center:
        try:
            chain_id, res_number, atom_name = user_center.split(":")
        except ValueError:
            raise cs.WrongAtomStringFormat(f"The specified atom is wrong '{user_center}'. \
Should be 'chain:resnumber:atomname'")
        for chain in structure.get_chains():
            if chain.id == chain_id:
                for residue in chain.get_residues():
                    if residue.id[1] == int(res_number):
                        for atom in residue.get_atoms():
                            if atom.name == atom_name: 
                                COI = np.array(list(atom.get_vector())) 
  
    # calculate protein and ligand COM
    com_protein = calculate_com(structure)
    com_ligand = calculate_com(ligand)

    # calculating the maximum d of the ligand
    coor_ligand = []
    for atom in ligand.get_atoms():
        coor_ligand.append(list(atom.get_vector() - com_ligand))

    coor_ligand = np.array(coor_ligand)
    coor_ligand_max = np.amax(coor_ligand, axis=0)
    d_ligand = np.sqrt(np.sum(coor_ligand_max ** 2))

    # set threshold for near and far contacts based on ligand d
    if d_ligand / 2 < 5.0:
        d5_ligand = 5.0
    else:
        d5_ligand = d_ligand / 2 + 1

    if d_ligand > 8.0:
        d8_ligand = d_ligand / 2 + 4
    else:
        d8_ligand = 8.0

    # calculate vector to move the ligandi
    if user_center:
        move_vector = com_ligand - COI
    else:
        move_vector = com_ligand - com_protein

    # translate the ligand to the protein COM (COI for PPI)
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
    D = 10.0 if user_center else np.ceil(6.0 + d)
    D_initial = D
    logger.info("Sampling {}A spherical box around the centre of the receptor/interface.".format(D))

    if user_center:
        sphere_cent = COI
    else:
        sphere_cent = com_protein

    j = 0
    logger.info("Generating {} poses...".format(nposes))
    start_time = time.time()
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
        new_ligand_COM = calculate_com(ligand)

        # rotate ligand
        vector = Vector(new_ligand_COM)
        rotation_matrix = rotaxis(np.random.randint(0, 2 * np.pi), vector)

        for atom in ligand.get_atoms():
            coords_after = atom.get_vector().left_multiply(rotation_matrix)
            atom.set_coord(coords_after)

        # check if it's inside the sampling sphere
        dist = np.sqrt((new_ligand_COM[0] - sphere_cent[0]) ** 2 + (new_ligand_COM[1] - sphere_cent[1]) ** 2 + (
                new_ligand_COM[2] - sphere_cent[2]) ** 2)

        if dist < D:
            # check contacts at: 5A (no contacts) and 8A (needs contacts)
            protein_list = Selection.unfold_entities(structure, "A")
            contacts5 = []
            contacts8 = []
            ligand_atoms = list(ligand.get_atoms())
            
            contacts5.append( NeighborSearch(protein_list).search(new_ligand_COM, d5_ligand, "S"))
            contacts8 = NeighborSearch(protein_list).search(new_ligand_COM, d8_ligand, "S")
            if contacts8 and not any(contacts5):
                j += 1
                io = PDBIO()
                io.set_structure(ligand)
                output_name = os.path.join(outputfolder, 'ligand{}.pdb'.format(j))
                io.save(output_name)
                output.append(output_name)
                start_time = time.time()

            end_time = time.time()
            total_time = end_time - start_time
            if total_time > 60:
                D += 1
                if D - D_initial >= 20:
                    logger.info("Original box increased by 20A. Aborting...")
                    break
                start_time = end_time
                logger.info("Increasing sampling box by 1A.")
    logger.info("{} poses created successfully.".format(j))
    return output, D, list(sphere_cent)


def join(receptor, ligands, residue, output_folder=".", output="input{}.pdb"):
    """
    Join receptor and ligand PDB files into one, conserving old formatting
    and not repeating atom numbers.
    """

    with open(receptor, "r") as f:
        lines = f.readlines()
        receptor_content = [line for line in lines if
                            line[17:20] != residue and not line.startswith("CONECT") and not line.startswith("END")]
        connects = [line for line in lines if line.startswith("CONECT")]
        ligand_content_without_coords = [line[0:27] + "{}" + line[56:] for line in lines if line[17:20] == residue]
        atom_nums = [line[6:11] for line in lines if line[17:20] == residue]

    outputs = []
    for i, ligand in enumerate(ligands):
        with open(ligand, "r") as fin:
            # exclude connects but keep initial atom names (CL problem)
            ligand_coords = {line[12:16].strip(): line[27:56] for line in fin if
                             line.startswith("ATOM") or line.startswith("HETATM")}
            # assert len(ligand_coords) == len(ligand_content_without_coords), "Experimental part - send an issue to github"

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
    parser.add_argument("--user_center", type=str, default=None,
                        help="Center of protein-protein interface (chain:residue number:atom name)")
    args = parser.parse_args()
    return os.path.abspath(args.ligand), os.path.abspath(
        args.receptor), args.resname, args.poses, args.output_folder, args.user_center


if __name__ == "__main__":
    ligand, receptor, resname, poses, output_folder, user_center = parse_args()
    output, sphere_cent, D = randomize_starting_position(ligand, receptor, output_folder, poses, None, user_center)
    # Make format ready to be captured"
    print("OUTPUT; {}. Sphere center: {}. D: {}".format(output, sphere_cent, D))
