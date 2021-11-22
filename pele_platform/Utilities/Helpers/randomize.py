from Bio.PDB import PDBParser, PDBIO, NeighborSearch, Selection, rotaxis
from Bio.PDB.vectors import Vector
import numpy as np
import os
import re
import time

import pele_platform.Errors.custom_errors as cs
from pele_platform.Utilities.Parameters.parameters import Parameters
from pele_platform.Utilities.Helpers import helpers


def calculate_com(structure):
    """
    Calculates center of mass of the structure (ligand or protein).

    Parameters
    ----------
    structure : biopython Structure object
        PDB of choice loaded into biopython (only chains of interest).

    Returns
    -------
        A list defining center of mass of the structure.
    """
    structure_mass = 0.0
    com = np.zeros(3)
    for atom in structure.get_atoms():
        com = com + np.array(list(atom.get_vector())) * atom.mass
        structure_mass += atom.mass
    com = com / structure_mass
    return com


def set_box(final_site, ligand_positions, system):
    """
    Sets the box radius based on the spread of initial poses generated during randomization and the final site
    coordinates.

    The following logic is applied:
        - the mid-point of all ligands is used as one end of the box (centroid)
        - user-define final site is used as the other end of the box (final_site)
        - box center is set as a mid-point between the centroid and the final site
        - the radius is the distance between the centroid and the final site plus 2 A.

    Parameters
    ----------
    final_site : str
        Final site set by the user, i.e. orthosteric_site (when running GPCR) or final_site (OutIn).
    ligand_positions : List[str]
        List of PDB files with initial ligand positions (before they are joined to the protein).
    system : str
        Path to PDB file with the protein-ligand complex.

    Returns
    -------
        Box radius, box_center.
    """
    parser = PDBParser()
    all_coms = []

    final_site = helpers.get_coords_from_residue(system, final_site)

    for file in ligand_positions:
        structure = parser.get_structure('input', file)
        all_coms.append(calculate_com(structure))

    centroid = np.average(all_coms, axis=0)  # middle point of where all the ligands were generated
    distance = np.linalg.norm(centroid - final_site)  # width of the whole box
    box_radius = distance / 2 + 2  # half of the with of the box + 2A
    box_center = np.average([final_site, centroid], axis=0)  # mid point between ligand centroid and final site
    return round(box_radius, 3), list(box_center)


def randomize_starting_position(parameters, box_center=None, box_radius=None):
    """
    Randomize initial ligand position around the receptor.

    Parameters
    ----------
    parameters : Parameters
        Parameters object passed from Adaptive.simulation.
    box_center : Union[List[float], str]
        Box center defined by the user, either using a list of coordinates or by specifying an atom, e.g. A:123:CA
    box_radius : float
        Box radius defined by the user.

    Returns
    -------
        A list of PDB files with random ligand positions.
    """
    np.random.seed(parameters.seed)

    # When running SiteFinder, the users can narrow down the area where inputs will be
    # spawned by setting box radius and center. We're setting box_center as COI to avoid crazy code changes until 2.0.
    if parameters.site_finder:
        parameters.center_of_interface = box_center if box_center else None

    # Retrieve atom information from the string (if PPI or SF)
    if parameters.center_of_interface:
        pattern = r"([A-z]\:\d{1,4}\:[A-Z0-9]{1,4})"
        match = re.match(pattern, parameters.center_of_interface)

        if not match:
            raise cs.WrongAtomStringFormat(f"The specified atom is wrong '{parameters.center_of_interface}'. "
                                           f"Should be 'chain:residue number:atom name'.")
        else:
            chain_id, res_number, atom_name = parameters.center_of_interface.split(":")

    # Load PDB files into biopython
    parser = PDBParser()
    structure = parser.get_structure('protein', parameters.receptor)
    ligand = parser.get_structure('ligand', parameters.ligand_ref)
    output = []
    COI = np.zeros(3)  # initializing to [0, 0, 0]

    # Get center of interface (if PPI or custom SF box)
    if parameters.center_of_interface:
        for chain in structure.get_chains():
            if chain.id == chain_id:
                for residue in chain.get_residues():
                    if residue.id[1] == int(res_number):
                        for atom in residue.get_atoms():
                            if atom.name == atom_name:
                                COI = np.array(list(atom.get_vector()))
                                print("Center of interface:", COI)

    # calculate protein and ligand COM
    com_protein = calculate_com(structure)
    com_ligand = calculate_com(ligand)

    # calculating the maximum dimensions of the ligand
    coor_ligand = []
    for atom in ligand.get_atoms():
        coor_ligand.append(list(atom.get_vector() - com_ligand))

    coor_ligand = np.array(coor_ligand)
    coor_ligand_max = np.amax(coor_ligand, axis=0)
    d_ligand = np.sqrt(np.sum(coor_ligand_max ** 2))

    # set threshold for near and far contacts based on ligand dimensions
    if d_ligand / 2 < 5.0:
        d5_ligand = 5.0
    else:
        d5_ligand = d_ligand / 2 + 1

    if d_ligand > 8.0:
        d8_ligand = d_ligand / 2 + 4
    else:
        d8_ligand = 8.0

    # calculate vector to move the ligand respectively to the center of interface of center of mass of the protein
    if parameters.center_of_interface:
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

    # Calculate radius of the sphere from the origin
    if parameters.center_of_interface:
        if parameters.site_finder and box_radius:
            D = box_radius
        elif box_radius:  # GPCR and OutIn
            D = box_radius
        else:  # ppi
            D = 10
    else:
        D = np.ceil(6.0 + d)

    D_initial = D
    parameters.logger.info("Sampling {}A spherical box around the centre of the receptor/interface.".format(D))

    if parameters.center_of_interface:
        sphere_cent = COI
    else:
        sphere_cent = com_protein

    j = 0
    parameters.logger.info("Generating {} poses...".format(parameters.poses))
    start_time = time.time()
    while j < parameters.poses:
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
            contacts5 = NeighborSearch(protein_list).search(new_ligand_COM, d5_ligand, "S")
            contacts8 = NeighborSearch(protein_list).search(new_ligand_COM, d8_ligand, "S")

            if contacts8 and not any(contacts5):
                j += 1
                io = PDBIO()
                io.set_structure(ligand)
                output_name = os.path.join(parameters.inputs_dir, 'ligand{}.pdb'.format(j))
                io.save(output_name)
                output.append(output_name)
                start_time = time.time()

            end_time = time.time()
            total_time = end_time - start_time
            if total_time > 60:
                D += 1
                if D - D_initial >= 20:
                    parameters.logger.info("Original box increased by 20A. Aborting...")
                    break
                start_time = end_time
                parameters.logger.info("Increasing sampling box by 1A.")
    parameters.logger.info("{} poses created successfully.".format(j))
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
