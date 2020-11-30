from dataclasses import dataclass
import glob
import numpy as np
import os
import pandas as pd

import pele_platform.Adaptive.simulation as si
from pele_platform.Analysis.plots import _extract_coords
from pele_platform.Utilities.Helpers import bestStructs as bs
from pele_platform.Utilities.Helpers.helpers import cd, parallelize
import pele_platform.Utilities.Parameters.pele_env as pv


@dataclass
class Simulation:
    """
    Base Simulation class to set parameters, define input for the next BuildingBlock, etc.
    Both input and output should always be EnviroBuilder type.

    One class to rule them all, one class to find them, one class to bring them all and in PELE bind them.
    """
    env: pv.EnviroBuilder

    def set_params(self, simulation_type):
        setattr(self.env, simulation_type, True)

    def finish_state(self, value):
        self.env.next_step = value

    def set_working_folder(self, folder_name):
        self.original_dir = os.path.abspath(os.getcwd())
        self.env.folder_name = folder_name
        if hasattr(self.env, "pele_dir"):
            self.env.initial_args.folder = os.path.join(os.path.dirname(self.env.pele_dir), self.env.folder_name)

    def create_folders(self):
        self.env.create_files_and_folders()


@dataclass
class GlobalExploration(Simulation):
    env: pv.EnviroBuilder
    folder_name: str

    def run(self):
        """
        Launch global exploration.
        """
        self.set_params(simulation_type="full")
        self.set_working_folder(self.folder_name)
        self.env.build_adaptive_variables(self.env.initial_args)
        self.create_folders()
        self.env = si.run_adaptive(self.env)
        self.finish_state(self.env.output)
        return self.env


#@dataclass
#class InducedFitFast(Simulation):
#    env: pv.EnviroBuilder
#    folder_name: str
#
#    def run(self):
#        """
#        Launch induced fit fast.
#        """
#        self.set_params(simulation_type="induced_fit_fast")
#        self.set_working_folder(self.folder_name)
#        simulation_params = si.run_adaptive(self.env)
#        self.finish_state(simulation_params.output)
#        return simulation_params


@dataclass
class InducedFitExhaustive(Simulation):
    """
    Launch induced fit exhaustive.
    """
    env: pv.EnviroBuilder
    folder_name: str

    def run(self):
        self.set_params(simulation_type="induced_fit_exhaustive")
        self.set_working_folder(self.folder_name)
        self.env.build_adaptive_variables(self.env.initial_args)
        self.create_folders()
        #self.env.system = self.env.next_step
        self.env.input = glob.glob(self.env.next_step)
        self.env = si.run_adaptive(self.env)
        self.finish_state(self.env.output)
        return self.env


#@dataclass
#class Rescoring(Simulation):
#    """
#    Launch rescoring simulation.
#    """
#    env: pv.EnviroBuilder
#    folder_name: str
#
#    def run(self):
#        self.set_params(simulation_type="rescoring")
#        self.set_working_folder(self.folder_name)
#        simulation_params = si.run_adaptive(self.env)
#        self.finish_state(simulation_params.output)
#        return simulation_params


@dataclass
class Selection:

    simulation_params: pv.EnviroBuilder
    folder: str
        
    def run(self):
        self.choose_refinement_input()
        return self.simulation_params

    def choose_refinement_input(self):
        """
        Choose input for refinement simulation.
        Scan top 75% binding energies, pick n best ones as long as ligand COMs are >= box radius away from each other.
        """
        simulation_path = os.path.join(self.simulation_params.pele_dir, self.simulation_params.output)
        n_best_poses = int(self.simulation_params.iterations * self.simulation_params.pele_steps * (
                self.simulation_params.cpus - 1) * 0.75)

        with cd(simulation_path):
            files_out, _, _, _, output_energy = bs.main(str(self.simulation_params.be_column), n_structs=n_best_poses,
                                                        path=".",
                                                        topology=self.simulation_params.topology,
                                                        logger=self.simulation_params.logger)

            snapshot = 0
            files_out = [os.path.join(self.simulation_params.pele_dir, "results", f) for f in files_out]
            input_pool = [[f, snapshot, self.simulation_params.residue, self.simulation_params.topology] for f in
                          files_out]
            all_coords = parallelize(_extract_coords, input_pool, self.simulation_params.cpus)
            coords = [list(c[0:3]) for c in all_coords]
            dataframe = pd.DataFrame(list(zip(files_out, output_energy, coords)),
                                     columns=["File", "Binding energy", "1st atom coordinates"])
            dataframe = dataframe.sort_values(["Binding energy"], ascending=True)

            inputs = self._check_ligand_distances(dataframe)
            directory = os.path.join(os.path.dirname(self.simulation_params.pele_dir), self.folder)   ######### this needs to be done with setting working folder like in Simulation I guess
            if not os.path.isdir(directory):
                os.makedirs(directory, exist_ok=True)
            for i in inputs:
                os.system("cp {} {}/.".format(i, directory))

            self.simulation_params.next_step = os.path.join(directory, "*.pdb")

    def _check_ligand_distances(self, dataframe):

        inputs = []
        input_coords = []
        distances = []
        n_inputs = self.simulation_params.cpus - 1

        while len(inputs) < n_inputs:
            for f, c in zip(dataframe['File'], dataframe['1st atom coordinates']):
                if not input_coords:
                    inputs.append(f)
                    input_coords.append(c)
                else:
                    for ic in input_coords:
                        distances.append(abs(np.linalg.norm(np.array(c) - np.array(ic))))
                    distances_bool = [d > self.simulation_params.box_radius for d in distances]
                    if all(distances_bool):
                        inputs.append(f)
                        input_coords.append(c)
            break  # make sure it stops after running out of files to check
        return inputs


@dataclass
class Pipeline:

    iterable: list
    env: pv.EnviroBuilder

    def run(self):
        output = []
        for simulation in self.iterable:
            folder = "{}_{}".format(self.iterable.index(simulation) + 1, simulation.__name__)
            self.env = simulation(self.env, folder).run()
            output.append(self.env)
        return output
