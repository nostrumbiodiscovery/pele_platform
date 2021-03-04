from dataclasses import dataclass
import os

from pele_platform.Utilities.Parameters import pele_env
from pele_platform.Utilities.Helpers import helpers
from pele_platform.Adaptive import simulation
from pele_platform.site_finder.cluster import cluster_best_structures


@dataclass
class CovalentDocking:
    env: pele_env.EnviroBuilder
    original_dir: str = os.path.abspath(os.getcwd())

    def run(self):
        """
        Runs the whole covalent docking pipeline.
        Returns
        -------
            A tuple of EnviroBuilder objects with job variables for both simulations.
        """
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

    def set_top_level_directory(self):
        """
        Sets top level working folder to contain all simulation steps.
        """
        working_folder = os.path.abspath("{}_Pele".format(self.env.residue))

        if not self.env.folder:
            self.working_folder = (
                helpers.is_repeated(working_folder)
                if not self.env.adaptive_restart
                else helpers.is_last(working_folder)
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

            with helpers.cd(output_path):
                cluster_best_structures(5, n_components=simulation1.n_components,
                                        residue=simulation1.residue, topology=simulation1.topology,
                                        directory=self.refinement_dir, logger=simulation1.logger)

    def set_refinement_perturbation_params(self):
        """
        Sets parameters for the refinement side chain perturbation, including the refinement distance (default = 10 A).
        """
        self.env.refinement_distance = self.env._refinement_distance if self.env._refinement_distance else 10.0
        self.env.folder = os.path.join(self.working_folder, "2_refinement")
        self.env.system = os.path.join(self.refinement_dir, "*.pdb")
