import os
import warnings
import subprocess

from pele_platform.constants import constants
from pele_platform.Errors import custom_errors


class PDBChecker:
    def __init__(self, file):
        """
        Initializes PDBChecker class and load relevant lines from the file.

        Parameters
        ----------
        file : str
            Path to PDB file.
        """
        self.file = file
        self.atom_lines, self.conect_lines = self._load_lines()

    def _load_lines(self):
        """
        Reads in the PDB file.

        Returns
        -------
        PDB lines.
        """
        with open(self.file, "r") as pdb_file:
            lines = pdb_file.readlines()

        atom_lines = [
            line
            for line in lines
            if line.startswith("ATOM") or line.startswith("HETATM")
        ]

        conect_lines = [line for line in lines if line.startswith("CONECT")]

        return atom_lines, conect_lines

    def check(self):
        """
        Performs all checks: protonation, negative residues and CONECT lines.

        Returns
        -------
        PDB file with added CONECT lines (necessary for Parametrizer).
        """
        self.check_protonation()
        self.check_negative_residues()
        self.file = self.check_conects()

        return self.file

    def check_protonation(self):
        """
        Checks for hydrogen atoms in the input PDB to ensure that it has been protonated.

        Raises
        ------
        ProtonationError if no hydrogen atoms are found in the input PDB.
        """
        hydrogen_lines = []

        for line in self.atom_lines:
            if line[12:16].strip().startswith("H"):
                hydrogen_lines.append(line)

        if len(hydrogen_lines) < 1:
            raise custom_errors.ProtonationError(
                "We did not find any hydrogen "
                "atoms in your system - looks "
                "like you forgot to "
                "protonate it."
            )

    def check_conects(self):
        """
        Checks the PDB file for CONECT lines and attempts to add them, if there are none.

        Returns
        --------
        PDB file with added CONECT lines.
        """
        if len(self.conect_lines) < 1:
            warnings.warn(
                f"PDB file {self.file} is missing the CONECT lines at the end!"
            )

            # Import and export with Schrodinger to add CONECT lines without making any other changes
            print("Adding CONECT lines with Schrodinger...")
            schrodinger_path = os.path.join(
                constants.SCHRODINGER, "utilities/prepwizard"
            )
            file_name, ext = os.path.splitext(self.file)
            conect_pdb_file = f"{file_name}_conect{ext}"
            command_pdb = f"{schrodinger_path} -nohtreat -noepik -noprotassign -noimpref -noccd -delwater_hbond_cutoff 0 -NOJOBID {self.file} {conect_pdb_file}"
            subprocess.call(command_pdb.split())

            return conect_pdb_file

        else:
            return self.file

    def check_negative_residues(self):
        """
        Checks if the PDB file contains negative residue numbers which can cause parsing errors.

        Raises
        -------
        IncorrectResidueNumbers
        """
        residue_numbers = [int(line[22:26].strip()) for line in self.atom_lines]

        if min(residue_numbers) < 0:
            raise custom_errors.IncorrectResidueNumbers(
                "PDB file contains negative residue numbers which are not supported. "
                "Please renumber them starting from 1."
            )
