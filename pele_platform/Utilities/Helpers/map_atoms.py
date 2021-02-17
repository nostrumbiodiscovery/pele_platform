from typing import Any, List

from pele_platform.Utilities.Helpers import helpers
from pele_platform.Utilities.Parameters import pele_env
from pele_platform.constants import constants


class AtomMapper:
    """
    If atom or residue numbers change during preprocessing, the AtomMapper will map them to the new PDB file based on
    atomic coordinates. All fields specified in constants.atom_string_flags are checked automatically.
    """

    def __init__(
        self,
        args: pele_env.EnviroBuilder,
        env: pele_env.EnviroBuilder,
        original_system: str,
    ) -> None:
        self.args = args
        self.logger = env.logger
        self.original_system = original_system
        self.preprocessed_pdb = env.system
        self.atom_string_flags = constants.atom_string_flags
        self.all_args = [
            arg
            for arg in self.atom_string_flags
            if getattr(self.args, arg, None) is not None
        ]

    def run(self) -> pele_env.EnviroBuilder:
        for arg in self.all_args:
            arg_value = getattr(self.args, arg)
            arg_value = atom_number_to_atom_string(self.original_system, arg_value)
            new_atom_string = self.check_atom_string(arg_value)
            setattr(self.args, arg, new_atom_string)
        return self.args

    def check_atom_string(self, args: List[str]) -> List[str]:
        """
        Checks if the atom string needs mapping by attempting to extract its coordinates.
        """
        output = []
        args = args if isinstance(args, list) else [args]
        for arg in args:
            try:
                helpers.get_coords_from_residue(self.original_system, arg)
                helpers.get_coords_from_residue(self.preprocessed_pdb, arg)
                output.append(arg)
            except Exception as e:
                self.logger.info("{} - mapping it now!".format(e))
                _, after = self.map_atom_string(
                    arg, self.original_system, self.preprocessed_pdb, self.logger
                )
                output.append(after)
        return output

    @staticmethod
    def map_atom_string(
        atom_string: str,
        original_input: str,
        preprocessed_file: str,
        logger=Any,
    ) -> (str, str):
        """
        Maps old atom string to a new atom string by comparing coordinates of the original and preprocessed PBD files.
        """

        # read in the original and preprocessed PDB lines
        with open(original_input, "r") as initial:
            initial_lines = initial.readlines()

        with open(preprocessed_file, "r") as prep:
            preprocessed_lines = prep.readlines()

        # retrieve atom info from the original PDB
        chain, resnum, atom_name = atom_string.split(":")

        # extract coordinates from the original PDB
        for i in initial_lines:
            if (
                (i.startswith("HETATM") or i.startswith("ATOM"))
                and i[21].strip() == chain
                and i[22:26].strip() == resnum
                and i[12:16].strip() == atom_name
            ):
                initial_coords = get_coords_from_line(i)

                # extract coordinates from preprocessed file and compare to the original one
                for p in preprocessed_lines:
                    preprocessed_coords = get_coords_from_line(p)

                    if preprocessed_coords == initial_coords:
                        new_atom_name, new_resnum, new_chain = get_atom_from_line(p)
                        before = "{}:{}:{}".format(chain, resnum, atom_name)
                        after = "{}:{}:{}".format(new_chain, new_resnum, new_atom_name)
                        logger.info("Atom {} mapped to {}.".format(before, after))
                        return before, after


def get_atom_from_line(line):
    """
    Extracts atom name, residue number and chain ID from a PDB line.
    """
    atom_name = line[12:16].strip()
    residue_number = line[22:26].strip()
    chain_id = line[21].strip()
    return atom_name, residue_number, chain_id


def get_coords_from_line(line):
    """
    Extracts atom coordinates from a PDB line based on chain ID, residue number and PDB atom name.
    """
    return line[30:54].split()


def atom_number_to_atom_string(pdb_file, number):
    """
    Converts PDB atom number to PELE's atom string format, if necessary.
    """
    if not isinstance(number, list):
        number = [number]

    output = []
    for n in number:
        if isinstance(n, int) or n.isdigit():
            with open(pdb_file, "r") as f:
                lines = f.readlines()
            for line in lines:
                if line[6:11].strip() == str(n) and (
                    line.startswith("HETATM") or line.startswith("ATOM")
                ):
                    atom_name, resnum, chain = get_atom_from_line(line)
                    output.append("{}:{}:{}".format(chain, resnum, atom_name))
        else:
            output.append(n)
    return output
