from collections import defaultdict
from prody.proteins.pdbfile import parsePDB
import pele_platform.constants.constants as constants
import pele_platform.Utilities.Helpers.map_atoms as map_atoms
from PPP.checks_module import CheckforGaps


terminal_constr_spring = 5
backbone_constr_spring = 0.5

CONSTR_ATOM = """{{ "type": "constrainAtomToPosition", "springConstant": {0}, "equilibriumDistance": 0.0, "constrainThisAtom": "{1}:{2}:{3}" }},"""
CONSTR_DIST = """{{ "type": "constrainAtomsDistance", "springConstant": {}, "equilibriumDistance": {}, "constrainThisAtom": "{}:{}:{}", "toThisOtherAtom": "{}:{}:{}" }},"""
CONSTR_CALPHA = """{{ "type": "constrainAtomToPosition", "springConstant": {2}, "equilibriumDistance": 0.0, "constrainThisAtom": "{0}:{1}:_CA_" }},"""


class AlphaConstraints(object):

    def __init__(self, pdb, interval, backbone_spring, terminal_spring):
        """
        Class to parse alpha carbons in the PDB file and constrain them according to the user-defined settings.

        Parameters
        ----------
        pdb: str
            path to the system PDB
        interval: int
            interval at which backbone CA are supposed to be constrainte
        backbone_spring: float
            spring constant for backbone CA constraints
        terminal_spring: float
            spring constant for terminal CA constraints
        """
        self.pdb = pdb
        self.interval = interval
        self.backbone_spring_constant = backbone_spring
        self.terminal_spring_constant = terminal_spring
        self.residues = self.get_all_residues()
        self.interval_residues = self._apply_interval()
        self.gaps = self.find_gaps()

    def build_constraints(self):
        """
        Run the whole carbon alpha constraints pipeline, including terminals, gaps and interval backbone.

        Returns:
            A list of string constraints ready to be injected into the PELE configuration file.
        """
        return self.constraints

    @staticmethod
    def add_constraints(chain, resnum, spring=0.5):
        """
        Allow the user to create their own constraints by passing chain, residue number and spring constant.

        Args:
            chain (str): chain ID, e.g. "B"
            resnum (str of int): residue number, e.g. 123
            spring (Union[int, float], default=0.5: spring constant to constrain the atom.

        Returns:
            A string to constrain selected carbon alpha with user-defined spring value.
        """
        return CONSTR_CALPHA.format(chain, resnum, spring)

    def get_all_residues(self):
        """
        Scans the protein and retrieves all alpha carbons.

        Returns:
            A dictionary of carbons where chain ID is the key and values are a list of residue numbers,
            e.g. {'A': [1, 2, 3]}.
        """
        output = defaultdict(list)

        with open(self.pdb, "r") as pdb_file:
            pdb_lines = pdb_file.readlines()

        pdb_lines = [line for line in pdb_lines if line.startswith("ATOM")]

        for line in pdb_lines:
            (
                atom_name,
                residue_number,
                residue_name,
                chain_id,
            ) = map_atoms.get_atom_from_line(line)

            # converting str to int to calculate intervals later
            residue_number = int(
                residue_number
            )
            if atom_name == "CA" and residue_name in constants.AMINO_ACIDS:
                output[chain_id].append(residue_number)

        return output

    def _apply_interval(self):
        """
        Takes a dict of all alpha carbons and returns only the ones that satisfy the interval requirement.

        Returns:
            An amended dictionary of backbone residues where chain ID is the key and values are a list of residue
            numbers
        """
        output = defaultdict(list)

        for chain, residue_numbers in self.residues.items():
            for number in residue_numbers[self.interval:-self.interval]:  # don't duplicate terminal constraints
                for already_present_id in output[chain]:
                    if abs(already_present_id - number) < self.interval:
                        break
                else:
                    output[chain].append(number)
        return output

    def find_gaps(self):
        """
        Scans all residues in the protein and checks whether the N of the current residue and the C of the previous one
        are within peptide bond distance (1.50 A) in order to find any gaps in the backbone.

        Returns:
            Two dictionaries, first one containing the residues involved in the gaps, second one with all the remaining
            residues. Each of them has the chain as key and the residues numbers (previous, current) involved in a bond
            as values.
        """
        structure = parsePDB(self.pdb)
        gaps, no_gaps = CheckforGaps(structure, 1.50)
        return gaps

    @property
    def constraints(self):
        """
        Fills out the constrain templates and makes sure the final output generates a valid JSON.

        Returns:
            A list of string constraints (backbone, gaps and terminal) ready to be injected into the PELE configuration
            file (pele.conf).
        """
        json_start = [
            """"constraints":[""",
        ]
        json_end = ["],"]

        all_constraints = (
                self.backbone_constraints
                + self.gaps_constraints
                + self.terminal_constraints
        )

        # strip the last comma to avoid invalid JSON
        all_constraints[-1] = all_constraints[-1].strip(
            ","
        )

        constraints = json_start + all_constraints + json_end
        return constraints

    @property
    def terminal_constraints(self):
        """
        Constrains the first and the last residue's carbon alpha in each chain with the default
        terminal_spring_constant.

        Returns:
            A list of carbon alpha (CA) constraints on the terminals.
        """
        output = []
        for chain, residue_numbers in self.residues.items():
            output.append(
                CONSTR_CALPHA.format(
                    chain, residue_numbers[0], self.terminal_spring_constant
                )
            )
            output.append(
                CONSTR_CALPHA.format(
                    chain, residue_numbers[-1], self.terminal_spring_constant
                )
            )

        return output

    @property
    def gaps_constraints(self):
        """
        Constraints terminal CAs of any gaps detected in the protein.

        Returns:
            A list of gap constraints.
        """
        gaps_constr = []
        for chain, residues in self.gaps.items():
            gaps_constr = [
                CONSTR_ATOM.format(terminal_constr_spring, chain, terminal, "_CA_")
                for terminals in residues
                for terminal in terminals
            ]
        return gaps_constr

    @property
    def backbone_constraints(self):
        """
        Constraints all alpha carbons extracted earlier.

        Returns:
            A list of backbone constraints at a specified interval (omitting the terminal CAs).
        """
        output = []

        for chain, residue_numbers in self.interval_residues.items():
            residue_numbers = list(set(residue_numbers))
            residue_numbers.sort()

            for residue in residue_numbers:
                output.append(
                    CONSTR_CALPHA.format(chain, residue, self.backbone_spring_constant)
                )

        return output


def retrieve_constraints(
        pdb_file,
        interval=10,
        back_constr=backbone_constr_spring,
        ter_constr=terminal_constr_spring,
):
    """
    Runs the full AlphaConstraints pipeline to return a list of constraints ready to inject into JSON.

    Parameters
    ----------
    pdb_file: str
        PDB file with the protein
    interval: int
        interval at which the backbone CAs are supposed to be constrained
    back_constr: float
        spring constant for the backbone constraints
    ter_constr: float
        spring constant for the terminal constraints

    Returns
    -------
        A list of string constraints ready to be injected into the PELE configuration file.
    """
    constr = AlphaConstraints(pdb_file, interval, back_constr, ter_constr)
    constraints = constr.build_constraints()
    return constraints
