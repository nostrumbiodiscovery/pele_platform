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


class ConstraintBuilder(object):
    def __init__(self, pdb, interval, backbone_spring, terminal_spring):
        self.pdb = pdb
        self.interval = interval
        self.backbone_spring_constant = backbone_spring
        self.terminal_spring_constant = terminal_spring
        self.residues = self.get_all_residues()
        self.interval_residues = self._apply_interval()
        self.gaps = self.find_gaps()

    def create_constraints(self):
        return self.constraints

    def get_all_residues(self):
        """
        Retrieves all alpha carbons in the protein.
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
        Takes a list of all alpha carbons and returns only the ones that satisfy the interval requirement.
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
        Finds gaps in the protein chain based on the distance between N of the current residue and the
        C of the previous one, uses PPP.
        """
        structure = parsePDB(self.pdb)
        output, _ = CheckforGaps(structure, 1.55)
        return output

    @property
    def constraints(self):
        """
        Fills out the constrain templates and makes sure the final output generates a valid JSON.
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

        output: List[str] - constraints lines ready to be injected into JSON
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
        Constraints all gaps.
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
    constr = ConstraintBuilder(pdb_file, interval, back_constr, ter_constr)
    constraints = constr.create_constraints()
    return constraints
