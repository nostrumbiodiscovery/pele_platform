from dataclasses import dataclass
import numpy as np
import os
import pandas as pd
import shutil

from pele_platform.Utilities.Helpers import bestStructs as bs
from pele_platform.Utilities.Helpers.helpers import cd, is_repeated, is_last, parallelize
from pele_platform.analysis.plot import _extract_coords
import pele_platform.Utilities.Parameters.parameters as pv
import pele_platform.Adaptive.simulation as si


@dataclass
class SiteFinderLauncher:
    args: pv.ParametersBuilder

    def run_site_finder(self) -> (pv.ParametersBuilder, pv.ParametersBuilder):
        """
        Launch site_finder simulation.
        1) Run global exploration to identify the most important pockets
        2) Run induced fit simulation to find deep pockets
        """

        self._set_params_global()
        self.global_simulation = self._launch_global()

        if not self.args.skip_refinement:
            self._choose_refinement_input()
            self._set_params_refinement()
            self._clean_temp_files()
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
            self.working_folder = (
                is_repeated(working_folder)
                if not self.args.adaptive_restart
                else is_last(working_folder)
            )
        else:
            self.working_folder = os.path.abspath(self.args.folder)
        self.args.folder = os.path.join(self.working_folder, "1_global_exploration")
        self.args.full = True  # needed for global exploration

    def _launch_global(self):
        sim_params = si.run_adaptive(self.args)
        return sim_params

    def _choose_refinement_input(self):
        """
        Scan top 75% best binding energies, pick n best ones as long as ligand COMs are >= box radius away from each other.
        """
        simulation_path = os.path.join(
            self.global_simulation.pele_dir, self.global_simulation.output
        )
        n_best_poses = int(
            self.global_simulation.iterations
            * self.global_simulation.pele_steps
            * (self.global_simulation.cpus - 1)
            * 0.25
        )

        if not self.args.debug:
            with cd(simulation_path):
                files_out, _, _, _, output_energy = bs.main(
                    str(self.args.be_column),
                    n_structs=n_best_poses,
                    path=".",
                    topology=self.global_simulation.topology,
                    logger=self.global_simulation.logger,
                )

                snapshot = 0
                files_out = [
                    os.path.join(
                        self.global_simulation.pele_dir,
                        self.global_simulation.output,
                        f,
                    )
                    for f in files_out
                ]
                input_pool = [
                    [
                        f,
                        snapshot,
                        self.global_simulation.residue,
                        self.global_simulation.topology,
                    ]
                    for f in files_out
                ]
                all_coords = parallelize(_extract_coords, input_pool, 1)
                coords = [list(c[0:3]) for c in all_coords]
                dataframe = pd.DataFrame(
                    list(zip(files_out, output_energy, coords)),
                    columns=["File", "Binding energy", "1st atom coordinates"],
                )
                self.dataframe = dataframe.sort_values(
                    ["Binding energy"], ascending=True
                )

                inputs = self._check_ligand_distances()
                directory = os.path.join(self.working_folder, "refinement_input")

                if not os.path.isdir(directory):
                    os.makedirs(directory, exist_ok=True)
                for i in inputs:
                    os.system("cp {} {}/.".format(i, directory))

    def get_n_inputs(self):
        maximum_inputs = self.global_simulation.cpus - 1
        n_inputs = 20 if maximum_inputs > 20 else maximum_inputs
        return n_inputs

    def _check_ligand_distances(self):

        inputs = []
        input_coords = []
        n_inputs = self.get_n_inputs()

        for file, coord in zip(
            self.dataframe["File"], self.dataframe["1st atom coordinates"]
        ):

            if (
                len(inputs) == n_inputs
            ):  # get out of the loop, if we have enough inputs already
                break

            if not input_coords:
                inputs.append(file)
                input_coords.append(coord)

            else:
                distances = []

                for ic in input_coords:
                    distances.append(
                        abs(np.linalg.norm(np.array(coord) - np.array(ic)))
                    )
                distances_bool = [d > 6 for d in distances]

                if all(distances_bool):
                    inputs.append(file)
                    input_coords.append(coord)

        return inputs

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
        self.args.box_center = self.global_simulation.box_center
        self.args.box_radius = self.global_simulation.box_radius

    def _clean_temp_files(self):
        output_dir = os.path.join(self.global_simulation.pele_dir, self.global_simulation.output)
        temp_bs_dir = os.path.join(output_dir, "BestStructs")

        if os.path.exists(temp_bs_dir):
            shutil.rmtree(temp_bs_dir)

    def _launch_refinement(self):

        with cd(self.original_dir):
            # with cd(self.working_folder):
            if not self.args.debug:
                sim_params = si.run_adaptive(self.args)
            else:
                sim_params = None

        return sim_params
