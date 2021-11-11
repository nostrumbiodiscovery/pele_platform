import os
import shutil
import warnings
import subprocess
import tempfile

from pele_platform.constants import constants
from pele_platform.Errors import custom_errors
from pele_platform.Utilities.Helpers import helpers


class PDBChecker:

    def __init__(self, file, pele_dir=None):
        """
        Initializes PDBChecker class and load relevant lines from the file.

        Parameters
        ----------
        file : str
            Path to PDB file.
        """
        self.file = file
        self.inputs_dir = os.path.join(pele_dir, "input") if pele_dir else os.getcwd()
        self.fixed_file = self.file
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
        Performs all checks: protonation, negative residues, capped termini and CONECT lines.
        Then copies the final output PDB to the current working directory with '_fixed' suffix.

        Returns
        -------
        PDB file with added CONECT lines (necessary for Parametrizer) and no capped termini.
        """
        with tempfile.TemporaryDirectory() as tmp_dir:
            with helpers.cd(tmp_dir):
                self.check_protonation()
                self.check_negative_residues()
                self.fixed_file = self.remove_capped_termini()
                added_conects_file = self.check_conects()

                if added_conects_file:
                    self.fixed_file = added_conects_file

                self.fixed_file = self.save_file()

        return self.fixed_file

    def save_file(self):
        """
        Copies the final file from the temporary directory to the inputs directory (e.g. LIG_Pele/input).

        Returns
        -------
            Path to the corrected file name.
        """
        if not os.path.exists(self.inputs_dir):
            os.mkdir(self.inputs_dir)

        fixed_filename = os.path.join(self.inputs_dir, os.path.basename(self.file).replace(".pdb", "_fixed.pdb"))
        shutil.copy(self.fixed_file, fixed_filename)

        return fixed_filename

    def remove_capped_termini(self):
        """
        Removes any lines containing ACE and NMA residues, then saves them to replace the original file.
        """
        temp_file = os.path.join(os.getcwd(), os.path.basename(self.file.replace(".pdb", "_nocaps.pdb")))

        to_remove = []
        for line in self.atom_lines:
            if line[17:20].strip() == "ACE" or line[17:20].strip() == "NMA":
                to_remove.append(line)

        if to_remove:
            warnings.warn(f"File {self.file} contains uncapped termini. Removing all ACE and NMA lines...")

        self.atom_lines = [line for line in self.atom_lines if line not in to_remove]

        with open(temp_file, "w+") as file:
            for line in self.atom_lines:
                file.write(line)
            for line in self.conect_lines:
                file.write(line)

        return temp_file

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
            conect_pdb_file = os.path.join(os.path.dirname(self.fixed_file),
                                           os.path.basename(self.fixed_file.replace(".pdb", "_conect.pdb")))
            command_pdb = f"{schrodinger_path} -nohtreat -noepik -noprotassign -noimpref -noccd -delwater_hbond_cutoff 0 -NOJOBID {self.fixed_file} {conect_pdb_file}"
            subprocess.call(command_pdb.split())

            return conect_pdb_file
        else:
            return None

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
