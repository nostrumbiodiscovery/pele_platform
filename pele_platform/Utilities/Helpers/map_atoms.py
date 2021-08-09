from typing import Any, List, Union

from pele_platform.Utilities.Helpers import helpers
from pele_platform.Utilities.Parameters import parameters
from pele_platform.constants import constants


class AtomMapper:
    """
    If atom or residue numbers change during preprocessing, the AtomMapper will map them to the new PDB file based on
    atomic coordinates. All fields specified in constants.atom_string_flags are checked automatically.

    Input
    args: pele_env.EnviroBuilder - initial arguments passed by the user, passed from Adaptive.simulation
    env: pele_env.EnviroBuilder - passed from Adaptive.simulation
    ppp_system: str - complex PDB file before any preprocessing (syst.system)
    flags_to_check: List[str] - list of YAML flags to check, default constants.atom_string_flags

    Output
    args: pele_env.EnviroBuilder - original arguments with overwritten atom strings wherever necessary
    """

    def __init__(
        self,
        args: parameters.ParametersBuilder,
        env: parameters.ParametersBuilder,
        original_system: str,
        flags_to_check: List[str] = None,
    ) -> None:
        self.args = args
        self.logger = env.logger
        self.original_system = original_system
        self.preprocessed_pdb = env.system
        self.atom_string_flags = (
            flags_to_check if flags_to_check else constants.atom_string_flags
        )
        self.all_args = [
            arg
            for arg in self.atom_string_flags
            if getattr(self.args, arg, None) is not None
        ]

    def run(self) -> parameters.ParametersBuilder:
        """
        Run the whole mapping process.

        Input
        self - AtomMapper instance

        Output
        args: pele_env.EnviroBuilder - the original user parameters with overwritten atom strings
        """
        for arg in self.all_args:
            arg_value = getattr(self.args, arg)
            arg_value = atom_number_to_atom_string(self.original_system, arg_value)
            new_atom_string = self.check_atom_string(arg_value)
            setattr(self.args, arg, new_atom_string)
        return self.args

    def check_atom_string(self, args: List[str]) -> List[str]:
        """
        Checks if the atom string needs mapping by attempting to extract its coordinates.

        Input
        args: List[str]

        Output
        output: List[str]
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
        logger: Any,
    ) -> (str, str):
        """
        Maps old atom string to a new atom string by comparing coordinates of the original and preprocessed PBD files.

        Parameters
        -----------
        atom_string : str
            Atom string following the 'chain:residue number:atom name' or residue string with 'chain:resnum' format
        original_input : str
            Path to PDB file before preprocessing (syst.system)
        preprocessed_file : str
            Path to PDB file after preprocessing (env.system)
        logger: Any

        Returns
        --------
        before, after: (str, str) - tuple containing old and new (mapped) atom string
        """

        # read in the original and preprocessed PDB lines
        with open(original_input, "r") as initial:
            initial_lines = [
                line
                for line in initial.readlines()
                if line.startswith("HETATM") or line.startswith("ATOM")
            ]

        with open(preprocessed_file, "r") as prep:
            preprocessed_lines = [line for line in prep.readlines() if line.startswith("HETATM") or line.startswith("ATOM")]

        # retrieve atom info from the original PDB
        try:
            chain, resnum, atom_name = atom_string.split(":")  # atom string
        except ValueError:
            chain, resnum = atom_string.split(":")  # residue string
            atom_name = None

        # extract coordinates from the original PDB
        initial_coords = None
        for initial_line in initial_lines:
            if (
                initial_line[21].strip() == chain
                and initial_line[22:26].strip() == resnum
                and (atom_name is None or initial_line[12:16].strip() == atom_name)
            ):
                initial_coords = get_coords_from_line(initial_line)

        # extract coordinates from preprocessed file and compare to the original one
        for p in preprocessed_lines:
            preprocessed_coords = get_coords_from_line(p)

            if initial_coords is not None and preprocessed_coords == initial_coords:
                new_atom_name, new_resnum, _, new_chain = get_atom_from_line(p)

                if atom_name is not None:
                    before = "{}:{}:{}".format(chain, resnum, atom_name)
                    after = "{}:{}:{}".format(new_chain, new_resnum, new_atom_name)
                else:
                    before = "{}:{}".format(chain, resnum)
                    after = "{}:{}".format(new_chain, new_resnum)

                logger.info("Atom {} mapped to {}.".format(before, after))
                return before, after


def get_atom_from_line(line: str) -> (str, str, str, str):
    """
    Extracts atom name, residue number and chain ID from a PDB line.

    Input
    line: str - PDB line

    Output
    atom_name: str - PDB atom name from the PDB line
    residue_number: str - residue number from the PDB line
    residue_name: str - residue name from the PDB line, e.g. SER
    chain_id: str - chain ID from the PDB line
    """
    atom_name = line[12:16].strip()
    residue_number = line[22:26].strip()
    residue_name = line[16:21].strip()
    chain_id = line[21].strip()
    return atom_name, residue_number, residue_name, chain_id


def get_coords_from_line(line):
    """
    Extracts atom coordinates from a PDB line based on chain ID, residue number and PDB atom name.

    Input
    line: str - PDB line

    Output
    string of coordinates extracted from the PDB line
    """
    return line[30:54].split()


def atom_number_to_atom_string(
    pdb_file: str, number: Union[int, List[int], List[str], str]
):
    """
    Converts PDB atom number to PELE's atom string format, if necessary.

    Input
    pdb_file: str - PDB file
    number: Union[int, List[int], List[str], str] - PDB atom numbers

    Output
    output: List[str] - list of atom strings following the 'chain:residue number:atom name' format.
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
                    atom_name, resnum, _, chain = get_atom_from_line(line)
                    output.append("{}:{}:{}".format(chain, resnum, atom_name))
        else:
            output.append(n)
    return output
