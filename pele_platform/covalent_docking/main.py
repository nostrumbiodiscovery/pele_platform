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
        self.job1 = simulation.run_adaptive(self.env)
        self.choose_refinement_input()
        self.set_refinement_perturbation_params()
        self.job2 = simulation.run_adaptive(self.env)
        return self.job1, self.job2

    def set_general_perturbation_params(self):
        """
        Sets parameters for the initial side chain perturbation, making sure we set the correct working folder and
        ignore refinement distance for now.
        """
        self.env.perturbation = False
        self.env._refinement_angle = self.env.refinement_angle
        self.env.refinement_angle = None
        self.set_top_level_directory()
        self.env.folder = os.path.join(self.working_folder, "1_covalent_docking")

        if isinstance(self.env.skip_ligand_prep, list):
            self.env.skip_ligand_prep.append(self.env.residue)
        else:
            self.env.skip_ligand_prep = self.env.residue

    def correct_system(self):
        """
        Moves the covalent ligand to the other residues, replaces HETATM with ATOM, then assigns the resulting PDB
        as the new system.
        """
        corrected_system = os.path.join(
            self.original_dir,
            os.path.basename(self.env.system.replace(".pdb", "_corrected.pdb")),
        )
        chain, residue_number = self.env.covalent_residue.split(":")

        pdb_corrector.run(
            self.env.system,
            chain,
            int(residue_number),
            corrected_system,
            ligand_resname=self.env.residue,
            ligand_chain=self.env.chain,
        )

        self.env.system = corrected_system

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

    def choose_refinement_input(self):
        """
        Extracts 1000 lowest binding energy structures and clusters them based on heavy atom ligand coordinates using
        Gaussian Mixture Model. A lowest energy representative from each cluster is selected as input for the refinement
        simulation.
        """
        self.refinement_dir = os.path.join(self.working_folder, "refinement_input")

        if not self.env.debug:
            output_path = os.path.join(self.job1.pele_dir, self.job1.output)

            analysis_object = analysis.Analysis(
                simulation_output=output_path,
                resname=self.job1.residue,
                chain=self.job1.chain,
                traj=self.job1.traj_name,
                topology=self.job1.topology,
                cpus=1,
                skip_initial_structures=False
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
        template_json = '"refinementAngle": {},'
        self.env.refinement_angle = (
            template_json.format(self.env._refinement_angle)
            if self.env._refinement_angle
            else template_json.format(10.0)
        )
        self.env.folder = os.path.join(self.working_folder, "2_refinement")
        self.env.system = os.path.join(self.refinement_dir, "cluster*.pdb")

    def get_residue_type(self):
        """
        Extracts name of the residue the covalent ligand is bound to before correcting the system.
        """
        chain, residue_number = self.env.covalent_residue.split(":")
        residue_type = helpers.get_residue_name(self.env.system, chain, residue_number)
        return residue_type
