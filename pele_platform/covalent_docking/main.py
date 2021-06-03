from dataclasses import dataclass
import glob
import os
import re

from pele_platform.Utilities.Parameters import parameters
from pele_platform.Utilities.Helpers import helpers
from pele_platform.Adaptive import simulation
from pele_platform.analysis import analysis

from frag_pele.Covalent import pdb_corrector


@dataclass
class CovalentDocking:
    env: parameters.Parameters
    original_dir: str = os.path.abspath(os.getcwd())
    refinement_dir: str = None
    working_folder: str = None
    job1: parameters.Parameters = None
    job2: parameters.Parameters = None

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

        if not self.env.debug:
            self.choose_refinement_input()
            self.set_refinement_perturbation_params()
            self.job2 = simulation.run_adaptive(self.env)
        else:
            self.job2 = None

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
            self.env.skip_ligand_prep = [self.env.residue]

    def correct_system(self):
        """
        Moves the covalent ligand to the other residues, replaces HETATM with ATOM, then assigns the resulting PDB
        as the new system. Then adds original CONECT lines back to the extracted LIG.pdb.
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

        self.retrieve_ligand_conects()
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
        n_inputs = int(self.job1.cpus / 6)
        max_top_clusters = n_inputs if n_inputs > 1 else 1  # tests only have 5 CPUs

        output_path = os.path.join(self.job1.pele_dir, self.job1.output)

        analysis_object = analysis.Analysis(
            simulation_output=output_path,
            resname=self.job1.residue,
            chain=self.job1.chain,
            traj=self.job1.traj_name,
            topology=self.job1.topology,
            cpus=1,
            skip_initial_structures=False,
        )

        analysis_object.generate_clusters(
            self.refinement_dir,
            clustering_type="meanshift",
            representatives_criterion="local_nonbonding_energy",
            max_top_clusters=max_top_clusters,
        )

    def set_refinement_perturbation_params(self):
        """
        Sets parameters for the refinement side chain perturbation, including the refinement distance (default = 10 A).
        """
        self.env.refinement_angle = self.env._refinement_angle
        self.env.folder = os.path.join(self.working_folder, "2_refinement")
        self.env.system = os.path.join(self.refinement_dir, "cluster*.pdb")
        self.env.no_ppp = True
        self.env.covalent_docking_refinement = True

        self.recover_templates_from_job1()

    def recover_templates_from_job1(self):
        """
        Sets templates created in the first part of the simulation as external templates and ensures
        the parametrization of those ligands is skipped during refinement.
        """
        templates = glob.glob(
            os.path.join(self.job1.pele_dir, "DataLocal/Templates/OPLS2005/Protein/*")
        ) + glob.glob(
            os.path.join(
                self.job1.pele_dir, "DataLocal/Templates/OPLS2005/HeteroAtoms/*"
            )
        )

        self.env.template = [
            template
            for template in templates
            if os.path.basename(template) != "templates_generated"
        ]

        self.env.rotamers = glob.glob(
            os.path.join(self.job1.pele_dir, "DataLocal/LigandRotamerLibs/*")
        )

        for template in self.env.rotamers:
            self.env.skip_ligand_prep.append(
                os.path.basename(template.strip(".rot.assign"))
            )

        self.env.skip_ligand_prep = list(set(self.env.skip_ligand_prep))

    def get_residue_type(self):
        """
        Extracts name of the residue the covalent ligand is bound to before correcting the system.
        """
        chain, residue_number = self.env.covalent_residue.split(":")
        residue_type = helpers.get_residue_name(self.env.system, chain, residue_number)
        return residue_type

    def retrieve_ligand_conects(self):
        """
        Maps atom numbers from the original system PDB and adds modified CONECT lines to the extracted ligand PDB.
        """
        original_atom_numbers = list()
        original_conects = list()

        covalent_chain, covalent_resnum = self.env.covalent_residue.split(":")
        extracted_ligand = os.path.join(os.getcwd(), f"{self.env.residue}.pdb")

        # Load original PDB
        with open(self.env.system, "r") as system_file:
            system_lines = [
                line
                for line in system_file.readlines()
                if line.startswith("ATOM")
                or line.startswith("HETATM")
                or line.startswith("CONECT")
            ]

        # Find ligand covalent residue lines in the original PDB and extract atom numbers
        for line in system_lines:
            if line[17:20].strip() == self.env.residue:
                original_atom_numbers.append(line[7:11].strip())
            if line[22:26].strip() == covalent_resnum and line[21] == covalent_chain:
                original_atom_numbers.append(line[7:11].strip())

        # Extract CONECT lines containing relevant atom numbers (residue and ligand)
        for line in system_lines:
            if line.startswith("CONECT") and any(
                number in line for number in original_atom_numbers
            ):
                original_conects.append(line)

        original_conects = "".join(original_conects)

        # Extract new ligand atom numbers
        with open(extracted_ligand, "r") as ligand_file:
            ligand_lines = ligand_file.readlines()
            new_atom_numbers = [line[7:11].strip() for line in ligand_lines]
            # Schrodinger needs 4 spaces, otherwise it makes a mess
            new_atom_numbers = [f"    {number}" for number in new_atom_numbers]

        # Go over CONECT lines and update atom numbers
        for old, new in zip(original_atom_numbers, new_atom_numbers):
            original_conects = re.sub(rf"\b{old}\b", new, original_conects)

        # Append mapped CONECT lines to the extracted LIG.pdb
        with open(extracted_ligand, "a") as new_ligand:
            new_ligand.write(original_conects)
