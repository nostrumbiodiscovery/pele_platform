from dataclasses import dataclass
import numpy as np
import os
import pandas as pd

import pele_platform.Adaptive.simulation as si
from pele_platform.Analysis.plots import _extract_coords
from pele_platform.Utilities.Helpers import bestStructs as bs
from pele_platform.Utilities.Helpers.helpers import is_repeated, is_last, cd, parallelize
import pele_platform.Utilities.Parameters.pele_env as pv


@dataclass
class GlobalExploration:

    args: pv.EnviroBuilder

    def run(self):
        """
        Launch global exploration.
        """
        self.set_params()
        self.set_working_folder()
        simulation_params = si.run_adaptive(self.args)
        return simulation_params

    def set_params(self):
        """
        Set parameters for global exploration.
        """
        self.args.full = True

    def set_working_folder(self):
        """
        Set working folder. Users can pick their own.
        """
        self.original_dir = os.path.abspath(os.getcwd())
        working_folder = os.path.abspath("{}_Pele".format(self.args.residue))
        if not self.args.folder:
            self.working_folder = is_repeated(working_folder) if not self.args.adaptive_restart else is_last(
                working_folder)
        else:
            self.working_folder = os.path.abspath(self.args.folder)
        self.args.folder = os.path.join(self.working_folder, "1_global_exploration")


@dataclass
class InducedFitFast:

    args: pv.EnviroBuilder

    def run(self):
        """
        Launch induced fit fast.
        """
        self.set_params()
        choose_refinement_input(self.args)
        simulation_params = si.run_adaptive(self.args)
        return simulation_params

    def set_params(self):
        """
        Set induced fit params making sure there are no arguments from the previous simulation (building block).
        """
        self.args.full = False
        self.args.induced_fit_fast = True


@dataclass
class InducedFitExhaustive:
    """
    Launch induced fit exhaustive.
    """
    args: pv.EnviroBuilder

    def run(self):
        self.set_params()
        self.set_working_folder()
        choose_refinement_input(self.args)
        simulation_params = si.run_adaptive(self.args)
        return simulation_params

    def set_params(self):
        self.args.full = None  # in case it was passed from previous building block
        self.args.poses = None
        self.args.induced_fit_exhaustive = True

        self.args.system = os.path.join(self.args.working_folder, "refinement_input/*.pdb")

        if not self.args.test:
            self.args.iterations = 20

    def set_working_folder(self):
        self.args.folder = os.path.join(self.args.working_folder, "2_induced_fit_exhaustive")


@dataclass
class Rescoring:

    args: pv.EnviroBuilder

    def run(self):
        self.set_params()
        choose_refinement_input(self.args)
        simulation_params = si.run_adaptive(self.args)
        return simulation_params

    def set_params(self):
        pass


def choose_refinement_input(simulation_params):
    """
    Choose input for refinement simulation.
    Scan top 75% binding energies, pick n best ones as long as ligand COMs are >= box radius away from each other.
    """
    simulation_path = os.path.join(simulation_params.pele_dir, simulation_params.output)
    n_best_poses = int(simulation_params.iterations * simulation_params.pele_steps * (
            simulation_params.cpus - 1) * 0.75)

    with cd(simulation_path):
        files_out, _, _, _, output_energy = bs.main(str(simulation_params.be_column), n_structs=n_best_poses, path=".",
                                                    topology=simulation_params.topology,
                                                    logger=simulation_params.logger)

        snapshot = 0
        files_out = [os.path.join(simulation_params.pele_dir, "results", f) for f in files_out]
        input_pool = [[f, snapshot, simulation_params.residue, simulation_params.topology] for f in
                      files_out]
        all_coords = parallelize(_extract_coords, input_pool, simulation_params.cpus)
        coords = [list(c[0:3]) for c in all_coords]
        dataframe = pd.DataFrame(list(zip(files_out, output_energy, coords)),
                                 columns=["File", "Binding energy", "1st atom coordinates"])
        dataframe = dataframe.sort_values(["Binding energy"], ascending=True)

        inputs = _check_ligand_distances(simulation_params, dataframe)
        directory = os.path.join(simulation_params.working_folder, "refinement_input")

        if not os.path.isdir(directory):
            os.makedirs(directory, exist_ok=True)
        for i in inputs:
            os.system("cp {} {}/.".format(i, directory))


def _check_ligand_distances(simulation_params, dataframe):

    inputs = []
    input_coords = []
    distances = []
    n_inputs = simulation_params.global_simulation.cpus - 1

    while len(inputs) < n_inputs:
        for f, c in zip(dataframe['File'], dataframe['1st atom coordinates']):
            if not input_coords:
                inputs.append(f)
                input_coords.append(c)
            else:
                for ic in input_coords:
                    distances.append(abs(np.linalg.norm(np.array(c) - np.array(ic))))
                distances_bool = [d > simulation_params.box_radius for d in distances]
                if all(distances_bool):
                    inputs.append(f)
                    input_coords.append(c)
        break  # make sure it stops after running out of files to check
    return inputs
