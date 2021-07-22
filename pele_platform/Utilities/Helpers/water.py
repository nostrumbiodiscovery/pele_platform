from dataclasses import dataclass, field
import glob
import numpy as np
import os
import re
from typing import List

from Bio.PDB import PDBParser, PDBIO, Selection, NeighborSearch, Vector

import pele_platform.constants.constants as cs
import pele_platform.Errors.custom_errors as ce
import pele_platform.Utilities.Helpers.helpers as hp
import pele_platform.constants.pele_params as pp

TEMPLATE = '{{"watersToPerturb": {{"links": {{"ids": [{index}] }}}}, "Box": {{"radius": {radius}, "fixedCenter": [{com}], "type": "sphericalBox"}}}}'


@dataclass
class WaterIncluder:
    input_pdbs: list
    n_waters: int
    user_waters: List = field(default_factory=lambda: [])
    ligand_perturbation_params: str = ""
    ligand_residue: str = ""
    water_center: bool = False
    water_radius: bool = False
    water_to_exclude: List = field(default_factory=lambda: [])
    sim_path: str = "."
    allow_empty_selectors: bool = False
    water_temp: int = 5000
    water_trials: int = 10000
    water_overlap: float = 0.78
    water_constr: float = 0
    water_freq: int = 1
    all_waters: str = ""
    perturbed_waters = []
    test: bool = False
    water_ids_to_track: List = field(default_factory=list)

    def run(self):
        """
        Runs the whole workflow to create water section in the control file.
        """
        self.set_empty_selectors()
        self.set_user_waters()
        self.water_checker()
        self.add_water()
        self.set_water_control_file()

    def set_empty_selectors(self):
        # If True include in control file else null
        self.allow_empty_selectors = (
            '"allowEmptyWaterSelectors": true,' if self.allow_empty_selectors else ""
        )

    def set_user_waters(self):
        """
        If the user set waters: "all_waters", it extracts them from the pdb.
        """
        if self.user_waters == "all_waters":
            self.user_waters = hp.retrieve_all_waters(self.input_pdbs[0])

    def set_parameters(self):
        """
        Adds parameters change JSON for waters.
        """
        if self.ligand_perturbation_params:
            self.ligand_perturbation_params = (
                self.ligand_perturbation_params.rstrip("]\n") + pp.WATER_PARAMS
            )
        else:
            self.ligand_perturbation_params = (
                ', "parametersChanges" : [ ' + pp.WATER_PARAMS.replace(",{", "{")
            )

    def set_box_center(self, inp=False):
        """
        Sets water box center around the ligand or user-defined water center (if no ligand).
        """
        if self.ligand_residue:
            com = ligand_com(inp, self.ligand_residue)
            com_format = ["{:.10f}".format(elem) for elem in com[0]]
            self.water_center = ", ".join(com_format)
        else:
            com_format = ["{:.10f}".format(coord) for coord in self.water_center]
            self.water_center = ", ".join(com_format)

    def set_box_radius(self):
        """
        Sets water box radius.
        """
        if self.water_radius:
            self.water_radius = self.water_radius
        else:
            self.water_radius = 6

    def set_water_input(self, inp):
        """
        Sets water sites JSON (box, radius, water IDs) by formatting the TEMPLATE.
        """
        try:
            self.set_box_center(inp)
        except TypeError:
            raise ce.NotCenterOfWaterBox(
                'Center of water box not found. Set it with "water_center" flag or ensure '
                "you specified the residue name, so it can be retrieved automatically."
            )
        self.set_box_radius()

        waters = sorted(hp.retrieve_all_waters(inp, exclude=self.water_to_exclude))
        self.water_ids_to_track = self.retrieve_indices_to_track(waters)
        waters = ", ".join(['"{}"'.format(water) for water in waters])
        self.all_waters += waters + ", "
        return TEMPLATE.format(
            index=waters, radius=self.water_radius, com=self.water_center
        )

    def set_water_control_file(self):
        """
        Formats the water perturbation string in constants.
        """
        if self.n_waters != 0 or self.user_waters:
            self.set_parameters()
            water_string = [
                self.set_water_input(inp).strip("'") for inp in self.input_pdbs
            ]
            self.all_waters = self.all_waters.rstrip(", ")
            self.water_line = cs.WATER.format(
                self.all_waters,
                self.allow_empty_selectors,
                self.water_temp,
                self.water_trials,
                self.water_overlap,
                self.water_constr,
                water_string,
            ).replace("'", "")
            self.water_energy = None  # TEMPORARY FIX

        else:
            self.water_energy = None
            self.water_radius = None
            self.water_center = None
            self.water_line = ""

    def water_checker(self):
        """
        Validates the number of water required to add by the user.

        Raises
        -------
        ValueError
            If number of water molecules is > 4 or < 1.
        """
        max_water = 4
        min_water = 1
        if self.n_waters:
            if self.n_waters > max_water:
                raise ValueError(
                    "Maximum {} water molecules are allowed.".format(max_water)
                )
            elif self.n_waters < min_water:
                raise ValueError(
                    "Number of water molecules (n_waters) has to be between {} and {}".format(
                        min_water, max_water
                    )
                )

    def add_water(self):
        """
        Adds water to the PDB file in a random position, then translates them until there are no clashes.
        """
        if self.test:
            np.random.seed(42)

        output = []
        n_inputs = len(self.input_pdbs)
        water_coords = []
        resnums = []
        atomnums = []
        chains = []
        resnames = []

        # Open the original PDB file
        with open(self.input_pdbs[0], "r") as file:
            # Figure out which lines refer to the actual structure and CONECTs, drop everything else
            lines = file.readlines()
            conect = [line for line in lines if "CONECT" in line]
            pdb_lines = [
                line for line in lines if "END" not in line and "CONECT" not in line
            ]

            for line in pdb_lines:
                if (
                    line.startswith("ATOM")
                    or line.startswith("HETATM")
                    or line.startswith("TER")
                ):
                    try:
                        # Extract atom information
                        resnum = line[22:27].strip()
                        atomnum = line[7:11].strip()
                        chain = line[21]
                        resname = line[17:20]
                        resnums.append(resnum)
                        atomnums.append(atomnum)
                        chains.append(chain)
                        resnames.append(resname)

                        # If there are already waters in the system but were not selected to be perturbed, we exclude
                        # them
                        if resname == "HOH":
                            water = f"{chain}:{resnum}"
                            if (
                                water not in self.user_waters
                                and water not in self.water_to_exclude
                            ):
                                self.water_to_exclude.append(water)
                    # Line too short - Remarks pdb
                    except IndexError:
                        pass

        # Return if no waters are supposed to be added
        if self.n_waters < 1:
            return

        else:
            # Check the maximum existing residue name, so we know where to introduce the waters
            lig_length = resnames.count(self.ligand_residue)
            resnums = [int(num) for num in resnums if num]
            max_resnum = max(resnums)
            water_resnums = []

            # Figure out the chain ID and atom numbers to introduce the waters
            water_chain = chains[0]  # water chain = 1st protein chain
            atomnum = max([int(num) for num in atomnums if num]) + 1 + lig_length

            # Enumerate enough water templates to add n_waters to each input
            water = cs.water * self.n_waters * n_inputs
            for input_pdb in range(n_inputs):
                for water_string in range(self.n_waters):
                    # Randomize oxygen coordinates - create an [x, y, z] vector
                    O_coords = Vector([np.random.randint(0, 100) for _ in range(3)])
                    # Add hydrogens to the oxygen
                    H1_coords = O_coords + Vector(0.757, 0.586, 0.0)
                    H2_coords = O_coords + Vector(-0.757, 0.586, 0.0)
                    water_coords = (
                        water_coords
                        + [list(O_coords)]
                        + [list(H1_coords)]
                        + [list(H2_coords)]
                    )
                    # Increment residue number, so each added water has a different one
                    max_resnum += 1
                    water_resnums = water_resnums + [max_resnum] * 3
                max_resnum += 1

            # Calculate atom numbers of all waters
            water_atomnums = [atomnum + j for j in range(self.n_waters * 3 * n_inputs)]

            # Create water PDB lines based on calculated atom numbers, residues, etc.
            water_output = []
            for atom, num, resnum, coord in zip(
                water, water_atomnums, water_resnums, water_coords
            ):
                # Format coordinates, so they fit into the PDB format
                coord = ["{:7.4f}".format(c) for c in coord]
                coord = " ".join(coord)
                water_output.append(atom.format(num, water_chain, resnum, coord))

            # Slice created water PDB lines and split between different input PDBs
            sliced_water_output = []
            for i in range(0, len(water_output), self.n_waters * 3):
                sliced_water_output.append(water_output[i : i + self.n_waters * 3])

            # Loop over PDB inputs and
            for input_pdb, water_output in zip(self.input_pdbs, sliced_water_output):
                new_protein_file = input_pdb
                # Write PDB lines followed by created water lines
                with open(input_pdb, "w+") as file:
                    for line in pdb_lines:
                        file.write(line)
                    file.write("\n")
                    for line in water_output:
                        file.write(line)
                    file.write("END")

                # Load the input PDB file again with Biopython to check for contacts
                parser = PDBParser()
                structure = parser.get_structure("complex", new_protein_file)
                water_list = []

                # Get all protein atoms to check for clashes
                protein_list = Selection.unfold_entities(structure, "A")

                # Get all relevant water atoms to check for clashes
                for res in structure.get_residues():
                    resnum = res._id[1]
                    if res.resname == "HOH":
                        if resnum not in resnums:
                            water_list = water_list + Selection.unfold_entities(
                                res, "A"
                            )

                # Check contacts between added waters and the protein at 5.0 angstrom
                contacts5 = []
                for water_output in water_list:
                    contacts5 = contacts5 + NeighborSearch(protein_list).search(
                        water_output.coord, 5.0, "A"
                    )
                contacts5 = [
                    c for c in contacts5 if c not in water_list
                ]  # exclude "self" contacts

                # Keep on tranlsating the water molecules as long as there are clashes at 5.0 A
                while contacts5:
                    contacts5 = []
                    for w_ in water_list:
                        x, y, z = w_.coord
                        # Set new coordinates and check contacts again
                        w_.set_coord([x - 5, y, z])
                        contacts5 = contacts5 + NeighborSearch(protein_list).search(
                            w_.coord, 5.0, "A"
                        )
                        contacts5 = [c for c in contacts5 if c not in water_list]

                # Save final output with translated water as a temporary file
                temp_protein_file = os.path.join(
                    os.path.dirname(input_pdb),
                    os.path.basename(input_pdb).replace(".pdb", "_temp.pdb"),
                )

                io = PDBIO()
                io.set_structure(structure)
                io.save(temp_protein_file)
                output.append(new_protein_file)

                # Open the temporary file created with biopython
                new_water_lines = []
                with open(temp_protein_file, "r") as temp:
                    temp_lines = temp.readlines()

                    # Iterate over lines created with biopython
                    for line in temp_lines:
                        if (
                            line[17:20].strip() == "HOH"
                            and int(line[22:27].strip()) not in resnums
                        ):
                            line = line.replace(
                                line[7:11], str(int(line[7:11]) + lig_length)
                            )
                            if line[12:15] == "2HW":
                                line = line.strip("\n") + "\nTER\n"
                            # If it's one of added waters, we manually change its residue number an save
                            new_water_lines.append(line)

                del new_water_lines[-1]  # Last biopython line is a not needed TER

                # Save new water lines, so they are not duplicated in the next run
                with open("added_waters.txt", "a+") as water_file:
                    for line in new_water_lines:
                        water_file.write(line)

                # Overwrite the original input PDB, save original PDB lines, added water lines and the original CONECTs.
                with open(new_protein_file, "w+") as file:
                    for line in pdb_lines:
                        file.write(line)
                    if not line.startswith("TER"):
                        file.write("TER\n")
                    for line in new_water_lines:
                        file.write(line)
                    for line in conect:
                        file.write(line)
                    file.write("\n")
                    file.write("END")

                # Remove temporary biopython file
                os.remove(temp_protein_file)

    @staticmethod
    def retrieve_indices_to_track(retrieved_waters):
        """
        Retrieves water indices to track in a manner compatible with Analysis.

        Parameters
        -----------
        retrieved_waters : List[str]
            List of user-defined or retrieved waters, e.g. ['A:202', 'A:203', 'A:204']

        Returns
        --------
        water_indices : List[tuple[str]]
            List of tuples indicating water indices, e.g. [("A", 202), ("A", 203), ("A", 204)]
        """
        output = list()

        for water in retrieved_waters:
            chain, resnum = water.split(":")
            output.append((chain, int(resnum)))

        return output


def water_ids_from_conf(configuration_file):
    """
    Extract IDs of perturbed water molecules from pele.conf file.

    Parameters
    -----------
    configuration_file : str
        Path to pele.conf file.

    Returns
    --------
    water_indices : List[tuple[str]]
        List of tuples indicating water indices, e.g. [("A", 202), ("A", 203), ("A", 204)]
    """
    with open(configuration_file, "r") as file:
        content = file.read()

    pattern = r"watersToPerturb\":.+?ids.+?(\[.+?\])"
    match = re.findall(pattern, content)

    if match:
        as_list = eval(match[0])
        output = WaterIncluder.retrieve_indices_to_track(as_list)
        # Easier to remove duplicates that do a super complicated regex
        return list(set(output))
    else:
        return []


def ligand_com(refinement_input, ligand_chain):
    """
    Calculate ligand's center of mass.

    Parameters
    ------------
    refinement_input : str
        Path to PDB file.

    ligand_chain : str
        Ligand chain ID.


    Returns
    --------
    output : list[float]
        Center of mass vector.
    """
    parser = PDBParser()
    output = []
    refinement_input = glob.glob(refinement_input)

    for inp in refinement_input:
        structure = parser.get_structure("inp", inp)
        mass = 0.0
        com = np.zeros(3)
        for res in structure.get_residues():
            if res.resname == ligand_chain:
                for atom in res.get_atoms():
                    com = com + np.array(list(atom.get_vector())) * atom.mass
                    mass += atom.mass
                com = com / mass

        output.append(com.tolist())

    return output
