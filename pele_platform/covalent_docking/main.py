from dataclasses import dataclass
import os
import subprocess

from pele_platform.Utilities.Parameters import parameters
from pele_platform.Utilities.Helpers import helpers
from pele_platform.Adaptive import simulation
from pele_platform.analysis import analysis
from pele_platform.constants import constants

from frag_pele.Covalent import pdb_corrector


@dataclass
class CovalentDocking:
    env: parameters.Parameters
    original_dir: str = os.path.abspath(os.getcwd())

    def run(self):
        """
        Runs the whole covalent docking pipeline.
        Returns
        -------
            A tuple of EnviroBuilder objects with job variables for both simulations.
        """
        self.env.residue_type = self.get_residue_type()
        self.correct_system()
        self.set_general_perturbation_params()
        job1 = simulation.run_adaptive(self.env)
        self.choose_refinement_input(job1)
        self.set_refinement_perturbation_params()
        job2 = simulation.run_adaptive(self.env)
        return job1, job2

    def set_general_perturbation_params(self):
        """
        Sets parameters for the initial side chain perturbation, making sure we set the correct working folder and
        ignore refinement distance for now.
        """
        self.env.perturbation = False
        self.env._refinement_distance = self.env.refinement_distance
        self.env.refinement_distance = None
        self.set_top_level_directory()
        self.env.folder = os.path.join(self.working_folder, "1_covalent_docking")

        if isinstance(self.env.skip_ligand_prep, list):
            self.env.skip_ligand_prep.append(self.env.residue)
        else:
            self.env.skip_ligand_prep = self.env.residue

    def correct_system(self):
        """
        Moves the covalent ligand to the other residues, replaces HETATM with ATOM and runs Protein Wizard to
        restore CONECT lines.
        """
        corrected_system = os.path.join(
            self.original_dir,
            os.path.basename(self.env.system.replace(".pdb", "_corrected.pdb")),
        )
        chain, residue_number = self.env.covalent_residue.split(":")
        schrodinger_output = corrected_system.replace("_corrected.pdb", "_final.pdb")

        pdb_corrector.run(
            self.env.system,
            chain,
            int(residue_number),
            corrected_system,
            ligand_resname=self.env.residue,
            ligand_chain=self.env.chain,
        )

        schrodinger_path = os.path.join(constants.SCHRODINGER, "utilities/prepwizard")
        command_pdb = "{} -nohtreat -noepik -noprotassign -noimpref -noccd -NOJOBID {} {}".format(schrodinger_path, corrected_system, schrodinger_output)
        os.system(command_pdb)

        if os.path.exists(schrodinger_output):
            self.env.system = schrodinger_output

    def set_top_level_directory(self):
        """
        Sets top level working folder to contain all simulation steps.
        """
        working_folder = os.path.abspath("{}_Pele".format(self.env.residue))

        if not self.env.folder:
            self.working_folder = (
                helpers.get_next_peledir(working_folder)
                if not self.env.adaptive_restart
                else helpers.get_next_peledir(working_folder)
            )
        else:
            self.working_folder = os.path.abspath(self.env.folder)

    def choose_refinement_input(self, simulation1):
        """
        Extracts 1000 lowest binding energy structures and clusters them based on heavy atom ligand coordinates using
        Gaussian Mixture Model. A lowest energy representative from each cluster is selected as input for the refinement
        simulation.

        Parameters
        ----------
        simulation1: pele_env.EnviroBuilder
            Job parameters of the initial simulation.
        """
        self.refinement_dir = os.path.join(self.working_folder, "refinement_input")

        if not self.env.debug:
            output_path = os.path.join(simulation1.pele_dir, simulation1.output)

            analysis_object = analysis.Analysis(
                simulation_output=output_path,
                resname=self.simulation1.residue,
                chain=self.simulation1.chain,
                traj=self.simulation1.traj_name,
                topology=self.simulation1.topology,
                cpus=1,
            )

            analysis_object.generate_clusters(
                self.refinement_dir,
                clustering_type="meanshift",
                representatives_criterion="local_nonbonding_energy",
            )

    def set_refinement_perturbation_params(self):
        """
        Sets parameters for the refinement side chain perturbation, including the refinement distance (default = 10 A).
        """
        template_json = '"refinementDistance": {}'
        self.env.refinement_distance = (
            template_json.format(self.env._refinement_distance)
            if self.env._refinement_distance
            else template_json.format(10.0)
        )
        self.env.folder = os.path.join(self.working_folder, "2_refinement")
        self.env.system = os.path.join(self.refinement_dir, "*.pdb")

    def get_residue_type(self):
        """
        Extracts name of the residue the covalent ligand is bound to before correcting the system.
        """
        chain, residue_number = self.env.covalent_residue.split(":")
        residue_type = helpers.get_residue_name(self.env.system, chain, residue_number)
        return residue_type
