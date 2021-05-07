from dataclasses import dataclass
import os
from pathlib import Path
import warnings

from pele_platform.Utilities.Helpers.helpers import (
    cd,
    get_next_peledir,
    get_latest_peledir,
)
import pele_platform.Utilities.Parameters.parameters as pv
import pele_platform.Adaptive.simulation as si


@dataclass
class SiteFinderLauncher:
    args: pv.ParametersBuilder

    def run_site_finder(self) -> (pv.ParametersBuilder, pv.ParametersBuilder):
        """
        Launch site_finder simulation.
        1) Run global exploration to identify the most important pockets.
        2) Run induced fit simulation to find deep pockets.

        Returns
        -------
        global_simulation: a ParametersBuilder object
            The ParametersBuilder object with global simulation parameters
        refinement_simulation: a ParametersBuilder object
            The ParametersBuilder object with local simulation parameters
        """
        self._set_params_global()
        self.global_simulation = self._launch_global()

        if not self.args.skip_refinement and not self.args.debug:
            self._generate_clusters()
            self._retrieve_external_templates()
            self._set_params_refinement()
            self.refinement_simulation = self._launch_refinement()
        else:
            self.refinement_simulation = None

        return self.global_simulation, self.refinement_simulation

    def _set_params_global(self):
        """
        Set parameters for global exploration. Users can choose their own
        working folder.
        """

        # Handling nested directories which will become pele_dir
        self.original_dir = os.path.abspath(os.getcwd())
        working_folder = os.path.abspath("{}_Pele".format(self.args.residue))

        if not self.args.folder:
            self.working_folder = (
                get_next_peledir(working_folder)
                if not self.args.adaptive_restart
                else get_latest_peledir(working_folder)
            )
        else:
            self.working_folder = os.path.abspath(self.args.folder)
        self.args.folder = os.path.join(self.working_folder,
                                        "1_global_exploration")

        # Global exploration parameters
        self.args.site_finder_global = True

        # Warn about using less than 60 CPUs
        if not self.args.test:
            if self.args.cpus < 60:
                warnings.warn(f"You are using only {self.args.cpus} "
                              f"CPUs. We recommend using at least 60 for "
                              f"the site_finder package.")

    def _launch_global(self):
        sim_params = si.run_adaptive(self.args)
        return sim_params

    def _generate_clusters(self):
        """
        Selects the system for refinement based on meanshift clustering and
        copies them into the refinement_input folder.
        """
        from pele_platform.analysis import analysis

        # Analysis parameters for refinement input selection
        n_inputs = int(self.global_simulation.cpus / 6)
        clustering_method = "meanshift"
        bandwidth = 10.0
        min_population = 0.001
        max_top_clusters = n_inputs if n_inputs > 1 else 1  # tests only have 5 CPUs
        top_clusters_criterion = "interaction_min"
        cluster_representatives_criterion = "interaction_min"

        # Directories within top-level "LIG_Pele"
        output_folder = os.path.join(self.global_simulation.folder,
                                     self.global_simulation.output)
        refinement_directory = os.path.join(self.working_folder,
                                            "refinement_input")

        # Instantiate Analysis and generate clusters
        analysis = analysis.Analysis(
            simulation_output=output_folder,
            resname=self.global_simulation.residue,
            chain=self.global_simulation.chain,
            traj=self.global_simulation.traj_name,
            topology=self.global_simulation.topology,
            cpus=1)

        analysis.generate_clusters(
            refinement_directory,
            clustering_type=clustering_method,
            bandwidth=bandwidth,
            min_population=min_population,
            max_top_clusters=max_top_clusters,
            top_clusters_criterion=top_clusters_criterion,
            representatives_criterion=cluster_representatives_criterion)

    def _set_params_refinement(self):
        """
        Sets parameters and working folder for the refinement simulation.
        """
        self.args.system = os.path.join(self.working_folder,
                                        "refinement_input/cluster*.pdb")
        self.args.folder = os.path.join(self.working_folder,
                                        "2_refinement_simulation")

        # Cleaning up args after global simulation
        self.args.site_finder_global = None
        self.args.poses = None

        # Refinement parameters
        self.args.site_finder_local = True
        self.args.box_center = self.global_simulation.box_center
        self.args.box_radius = self.global_simulation.box_radius

    def _retrieve_external_templates(self):
        """
        Finds templates and rotamers created in the global exploration and adds them to the refinement input arguments,
        then ensures they will not be parametrized again by adding them to "skip_ligand_prep" flag.
        """
        datalocal_dir = os.path.join(self.global_simulation.pele_dir, "DataLocal")
        templates = [str(path) for path in Path(datalocal_dir).rglob('*z')]
        rotamers = [str(path) for path in Path(datalocal_dir).rglob('*.assign')]

        self.args.template = templates
        self.args.rotamers = rotamers

        templates_names = [os.path.basename(file).rstrip("z") for file in templates]
        rotamers_names = [os.path.basename(file).rstrip("rot.assign") for file in rotamers]

        self.args.skip_ligand_prep = list(set(templates_names + rotamers_names))

    def _launch_refinement(self):
        """
        Launches refinement simulation (local), if not in debug mode.
        """
        with cd(self.original_dir):
            if not self.args.debug:
                sim_params = si.run_adaptive(self.args)
            else:
                sim_params = None

        return sim_params
