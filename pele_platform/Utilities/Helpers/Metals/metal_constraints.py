import pele_platform.constants.constants as cs
import pele_platform.Errors.custom_errors as ce
from pele_platform.Utilities.Helpers import map_atoms
from Bio.PDB import PDBParser, NeighborSearch, Selection, Vector, vectors
import itertools
import numpy as np


def find_metals(protein_file):

    # read in the protein file
    parser = PDBParser()
    structure = parser.get_structure("protein", protein_file)

    # find metals
    metals = []
    for chain in structure.get_chains():
        for residue in chain.get_residues():
            for atom in residue.get_atoms():
                if atom.element in cs.metals:
                    metals.append([atom, residue, chain])
    return metals, structure


def map_metal_constraints(protein_file, original_input, original_constraints, logger):

    atoms = []

    for orig in original_constraints:
        try:
            k, dist, atom1, atom2 = orig.split("-")
        except ValueError:  # If more than one atom
            continue
        atoms.extend([atom1, atom2])

    for atom in atoms:
        before, after = map_atoms.AtomMapper.map_atom_string(
            atom, original_input, protein_file, logger=logger
        )

        for j in range(len(original_constraints)):
            for b, a in zip(before, after):
                original_constraints[j] = original_constraints[j].replace(b, a)

    return original_constraints


def angle_classification(combinations, permissive):

    # angle classification
    ang_180 = []
    ang_90 = []
    ang_109 = []
    coordinated_atoms = []

    if permissive:
        lower = 0.65
        upper = 1.35
    else:
        lower = 0.8
        upper = 1.2

    for c in combinations:
        a = c[2]
        if 180 * lower <= a <= 180 * upper:
            ang_180.append(c)
        if 90 * lower <= a <= 90 * upper:
            ang_90.append(c)
        if 109.5 * lower <= a <= 109.5 * upper:
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

    return geo, coordinated_atoms


def find_geometry(
    metals, structure, permissive=False, all_metals=False, external=None, logger=None
):

    # check metal contacts
    output = []
    checked_metals = []
    structure_list = Selection.unfold_entities(structure, "A")

    for metal in metals:

        # search distance based on metal type
        if metal[0].element == "YB":
            dist = 3.5
        elif metal[0].element == "K":
            dist = 3.3
        else:
            dist = 2.9

        metal_str = "{}:{}:{}".format(metal[2].id, metal[1].get_id()[1], metal[0].name)
        in_ext = []

        for i in external:
            if metal_str in i:
                in_ext = True

        if not in_ext and list(metal[0].coord) not in checked_metals:
            coords = metal[0].coord

            contacts = []

            for chain in structure.get_chains():

                for residue in chain.get_residues():
                    contacts_atoms = NeighborSearch(structure_list).search(
                        coords, dist, "A"
                    )
                    # exclude self-contacts, carbons and hydrogens
                    excluded_contacts = cs.metals + ["C", "H"]
                    contacts_atoms = [
                        c for c in contacts_atoms if c.element not in excluded_contacts
                    ]

                    for atom in contacts_atoms:
                        if (
                            residue in chain.get_residues()
                            and atom in residue.get_atoms()
                        ):
                            contacts.append([atom, residue, chain])

            combinations = list(itertools.combinations(contacts, 2))
            combinations = [list(c) for c in combinations]

            # get all atom - metal - atom angles
            for c in combinations:
                vi = Vector(c[0][0].coord)
                vj = Vector(c[1][0].coord)
                angle = vectors.calc_angle(vi, coords, vj) * 180 / np.pi
                c.append(angle)

            geo, coordinated_atoms = angle_classification(combinations, False)

            if geo is None and permissive:
                geo, coordinated_atoms = angle_classification(combinations, True)

                if geo is None and all_metals and combinations:

                    geo, coordinated_atoms = angle_classification(combinations, True)
                    if geo:
                        logger.info(
                            "Found {} geometry around {} (residue {}). Adding constraints.".format(
                                geo, metal[0].name, metal[1].get_id()[1]
                            )
                        )
                        checked_metals.append(list(metal[0].coord))
                    else:
                        coordinated_atoms = combinations
                        checked_metals.append(list(metal[0].coord))
                        geo = "no"
                        logger.info(
                            "Found {} geometry around {} (residue {}). Adding constraints to all atoms within {}A of the metal.".format(
                                geo, metal[0].name, metal[1].get_id()[1], dist
                            )
                        )

                elif geo is None and not all_metals:
                    raise ce.NoGeometryAroundMetal(
                        "Failed to determine geometry around {} (residue {}). Add constraints manually or set 'constrain_all_metals: true' to constrain all atoms within {}A of the metal.".format(
                            metal[0].name, metal[1].get_id()[1], dist
                        )
                    )

                elif geo is None and all_metals and not combinations:
                    logger.info(
                        "No atoms coordinated to {} (residue {}).".format(
                            metal[0].name, metal[1].get_id()[1]
                        )
                    )

                elif geo:
                    checked_metals.append(list(metal[0].coord))
                    logger.info(
                        "Found {} geometry around {} (residue {}). Adding constraints.".format(
                            geo, metal[0].name, metal[1].get_id()[1]
                        )
                    )

            elif geo is None and all_metals and combinations:
                geo, coordinated_atoms = angle_classification(combinations, True)

                if geo is None:
                    geo = "no"
                    coordinated_atoms = combinations
                    checked_metals.append(list(metal[0].coord))
                    logger.info(
                        "Found {} geometry around {} (residue {}). Adding constraints to all atoms within {}A of the metal.".format(
                            geo, metal[0].name, metal[1].get_id()[1], dist
                        )
                    )

                else:
                    logger.info(
                        "Found {} geometry around {} (residue {}). Adding constraints.".format(
                            geo, metal[0].name, metal[1].get_id()[1]
                        )
                    )

            elif geo is None and all_metals and not combinations:
                logger.info(
                    "No atoms coordinated to {} (residue {}).".format(
                        metal[0].name, metal[1].get_id()[1]
                    )
                )

            elif geo is None and not all_metals and not permissive:

                if not combinations:
                    logger.info(
                        "No atoms coordinated to {} (residue {}).".format(
                            metal[0].name, metal[1].get_id()[1]
                        )
                    )
                else:
                    geo, coordinated_atoms = angle_classification(combinations, True)

                    if geo is None:
                        geo = "no"
                        coordinated_atoms = combinations
                        checked_metals.append(list(metal[0].coord))
                        logger.info(
                            "Found {} geometry around {} (residue {}). Adding constraints to all atoms within {}A of the metal.".format(
                                geo, metal[0].name, metal[1].get_id()[1], dist
                            )
                        )
            else:
                checked_metals.append(list(metal[0].coord))
                logger.info(
                    "Found {} geometry around {} (residue {}). Adding constraints.".format(
                        geo, metal[0].name, metal[1].get_id()[1]
                    )
                )

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
                out = yaml_string.format(
                    spring_const,
                    atom_dist,
                    chain1,
                    resnum1,
                    atomname1,
                    chain2,
                    resnum2,
                    atomname2,
                )

                output.append(out)

            output = list(set(output))

            if output:
                output = ["{}".format(o) for o in output]

    return output


def main(
    original_constraints,
    protein_file,
    original_input,
    permissive=False,
    all_metals=False,
    external=None,
    logger=None,
):
    metals, structure = find_metals(protein_file)
    if external:
        external = map_metal_constraints(
            protein_file, original_input, original_constraints, logger=logger
        )
    output = find_geometry(
        metals, structure, permissive, all_metals, external, logger=logger
    )
    return output, external
