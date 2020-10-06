from dataclasses import dataclass
import os
from pele_platform.Allosteric.cluster import cluster_best_structures
from pele_platform.Utilities.Helpers.helpers import cd, is_repited, is_last
import pele_platform.Utilities.Parameters.pele_env as pv
import pele_platform.Adaptive.simulation as si


@dataclass
class AllostericLauncher:
    args: pv.EnviroBuilder

    def run_allosteric(self) -> (pv.EnviroBuilder, pv.EnviroBuilder):
        """
        Launch allosteric simulation.
        1) Run global exploration to identify the most important pockets
        2) Run induced fit simulation to find deep pockets
        """

        self._set_params_global()
        self.global_simulation = self._launch_global()

        if not self.args.skip_refinement:
            self._choose_refinement_input()
            self._set_params_refinement()
            self.refinement_simulation = self._launch_refinement()
        else:
            self.refinement_simulation = None

        return self.global_simulation, self.refinement_simulation

    def _set_params_global(self):
        """
        Set parameters for global exploration. Users can choose their own working folder.
        """
        self.original_dir = os.path.abspath(os.getcwd())
        working_folder = os.path.abspath("{}_Pele".format(self.args.residue))
        if not self.args.folder:
            self.working_folder = is_repited(working_folder) if not self.args.adaptive_restart else is_last(
                working_folder)
        else:
            self.working_folder = os.path.abspath(self.args.folder)
        self.args.folder = os.path.join(working_folder, "1_global_exploration")

        self.args.full = True  # needed for global exploration

    def _launch_global(self):
        sim_params = si.run_adaptive(self.args)
        return sim_params

    def _choose_refinement_input(self):
        simulation_path = os.path.join(self.global_simulation.pele_dir, self.global_simulation.output)

        if not self.args.debug:
            with cd(simulation_path):
                cluster_best_structures("5", n_components=self.global_simulation.n_components,
                                        residue=self.global_simulation.residue,
                                        topology=self.global_simulation.topology,
                                        directory=self.working_folder, logger=self.global_simulation.logger)

    def _set_params_refinement(self):

        self.args.system = os.path.join(self.working_folder, "refinement_input/*.pdb")
        self.args.folder = os.path.join(self.working_folder, "2_refinement_simulation")
        self.args.full = None
        self.args.poses = None
        self.args.induced_fit_exhaustive = True
        self.args.box_center = self.global_simulation.box_center
        self.args.box_radius = self.global_simulation.box_radius

        if not self.args.test:
            self.args.iterations = 20
            self.args.pele_steps = 10
        self.args.box_center = self.global_simulation.box_center
        self.args.box_radius = self.global_simulation.box_radius

    def _launch_refinement(self):

        with cd(self.original_dir):
            if not self.args.debug:
                sim_params = si.run_adaptive(self.args)
            else:
                sim_params = None

        return sim_params
