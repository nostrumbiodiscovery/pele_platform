import glob
import os
import subprocess
from typing import List

import PPP.main as ppp

from pele_platform.Adaptive.simulation import run_adaptive
from pele_platform.Utilities.Parameters.parameters import Parameters
from pele_platform.constants import constants
from pele_platform.Utilities.Helpers import helpers


MANDATORY_FLAGS = ["chain", "data", "input"]

RESIDUE_MAPPING = {
    "HSD": "HID",
    "HSE": "HIE",
}


class ProteinProteinRefinement:
    """
    Protein-protein refinement workflow aims to optimize the interactions between two proteins using a modified induced
    fit protocol. It should manage inputs from: MaSIF, ClusPro, pydock, zdock, LightDock and Rosetta.
    """

    def __init__(self, env: Parameters):
        """
        Initializes ProteinProteinRefinement class.

        Parameters
        ----------
        env : Parameters
            Parameters object with all user-defined arguments and defaults.
        """
        self.env = env
        self.env.folder = self.set_working_folder()
        if not self.env.input:
            self.env.input = glob.glob(self.env.system)
        self.temp_files = []

    def run(self) -> Parameters:
        """
        Runs the whole protein-protein docking workflow.

        Returns
        -------
            Parameters object containing all job settings.
        """
        self.check_mandatory_flags()
        self.system_preprocessing()
        job_parameters = run_adaptive(self.env)
        self.clean_up()
        return job_parameters

    def set_working_folder(self) -> str:
        """
        Sets the top level directory in a regular way, but passing "ProteinProtein_Pele" directly, instead of generating
        the directory name based on the resname.

        Returns
        -------
            String with a non-repeated PELE directory.
        """

        main_dir = "ProteinProtein_Pele"

        if self.env.restart or self.env.adaptive_restart or self.env.only_analysis:
            return helpers.get_latest_peledir(main_dir)
        else:
            return helpers.get_next_peledir(main_dir)

    def check_mandatory_flags(self) -> None:
        """
        Ensures the user has provided all the necessary arguments.
        """
        errors = []

        for flag in MANDATORY_FLAGS:
            if not hasattr(self.env, flag):
                errors.append(flag)

        if errors:
            raise ValueError(f"The input file is missing some arguments: {errors}.")

    def system_preprocessing(self) -> None:
        """
        Fixes PDBs, runs Schrodinger and PPP.
        """
        self.env.input = self.pdb_cleanup()
        self.env.input = self.schrodinger_preprocessing()
        self.env.input = self.ppp_preprocessing()

        # We need to update the system, in case PPP renumbered residues, in which case all parameters generated based
        # on system would fail. This is not how stuff should be working though...
        # TODO: Improve this for 2.0 when we get rid of Adaptive.simulation nonsense.
        self.env.system = self.env.input[0]

    def pdb_cleanup(self) -> List[str]:
        """
        Removes any artifacts introduced to the PDB file by another software, e.g. HEADER and END lines in the middle
        of the protein.

        Returns
        -------
            A list of preprocessed PDB files.
        """
        output_list = []

        for input_file in self.env.input:

            output_file = os.path.join(
                os.getcwd(), os.path.basename(input_file).replace(".pdb", "_fix.pdb")
            )

            with open(input_file, "r") as f:

                # Remove HEADER, END lines and any other random stuff
                lines = [
                    line
                    for line in f.readlines()
                    if line.startswith("ATOM") or line.startswith("HETATM")
                ]

                for idx, line in enumerate(lines):

                    # Convert His to Schrodinger nomenclature
                    for key, value in RESIDUE_MAPPING.items():
                        lines[idx] = lines[idx].replace(key, value)
                        lines[idx] = lines[idx].replace(key, value)

                    # Remove RA0 strings from cluspro based on the index (otherwise prody cannot find the atom type)
                    lines[idx] = lines[idx][0:68] + "\n"

            with open(output_file, "w") as fout:
                print(f"Cleaning up: {input_file}.")
                for line in lines:
                    fout.write(line)

            output_list.append(os.path.abspath(output_file))

        self.temp_files.extend(output_list)
        return output_list

    def schrodinger_preprocessing(self) -> List[str]:
        """
        Runs Schrodinger Protein Preparation Wizard with the following settings:
            - do not fill loops or side chains
            - do not run PROPKA
            - do not run protassign to optimize hydrogen bonding
            - delete waters
            - run without job control.

        Returns
        -------
            A list of preprocessed PDB files.
        """
        schrodinger_path = os.path.join(
            constants.SCHRODINGER, "utilities", "prepwizard"
        )
        schrodinger_command = (
            "{} -delwater_hbond_cutoff 5 -nopropka -noprotassign {} {} -NOJOBID"
        )

        output_list = []

        for input_file in self.env.input:
            print(f"Running Schrodinger Protein Preparation Wizard: {input_file}.")
            output = os.path.basename(input_file).replace(".pdb", "_sch.pdb")
            command = schrodinger_command.format(schrodinger_path, input_file, output)
            subprocess.call(command.split())
            output_list.append(output)

        self.temp_files.extend(output_list)
        return output_list

    def ppp_preprocessing(self) -> List[str]:
        """
        Runs PPP.

        Returns
        -------
            A list of preprocessed PDB files.
        """
        output_list = []

        for input_file in self.env.input:
            print(f"Running PPP: {input_file}.")
            preprocessed_input, missing_residues, _, _, _ = ppp.main(
                input_file,
                ".",
                charge_terminals=self.env.charge_ter,
                no_gaps_ter=self.env.gaps_ter,
            )
            output_list.append(preprocessed_input)

        return output_list

    def clean_up(self) -> None:
        """
        Cleans up all the temporary "mae", "sch" and "fix" files generated during preprocessing.
        """
        print("Cleaning up temporary files.")

        epik_files = glob.glob("*epik*")
        impref_files = glob.glob("*impref*")

        for file in self.temp_files + epik_files + impref_files:
            if os.path.isfile(file):
                os.remove(file)
