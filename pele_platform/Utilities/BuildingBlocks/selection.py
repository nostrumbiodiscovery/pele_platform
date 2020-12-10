from dataclasses import dataclass
from multiprocessing import Pool
import numpy as np
import os
import pandas as pd
from sklearn.mixture import GaussianMixture

from pele_platform.Analysis.plots import _extract_coords
from pele_platform.Utilities.Helpers import bestStructs as bs
from pele_platform.Utilities.Helpers.helpers import cd, parallelize
from pele_platform.Utilities.Parameters import pele_env as pv


@dataclass
class Selection:
    """
    Base class to handle all input selection algorithms, copy files, set next_step, etc.
    """
    simulation_params: pv.EnviroBuilder

    def copy_files(self):
        if not os.path.isdir(self.directory):
            os.makedirs(self.directory, exist_ok=True)
        for i in self.inputs:
            os.system("cp {} {}/.".format(i, self.directory))

    def set_next_step(self):
        self.simulation_params.next_step = os.path.join(self.directory, "*.pdb")

    def extract_poses(self, percentage, n_poses=None):
        simulation_path = os.path.join(self.simulation_params.pele_dir, self.simulation_params.output)
        n_best_poses = n_poses if n_poses else int(
            self.simulation_params.iterations * self.simulation_params.pele_steps * (
                    self.simulation_params.cpus - 1) * percentage)
        with cd(simulation_path):
            files_out, _, _, _, output_energy = bs.main(str(self.simulation_params.be_column), n_structs=n_best_poses,
                                                        path=".", topology=self.simulation_params.topology,
                                                        logger=self.simulation_params.logger)
        
        files_out = [os.path.join(self.simulation_params.pele_dir, "results", f) for f in files_out]
        return files_out, output_energy

    def extract_all_coords(self, files_out):
        snapshot = 0
        pool = Pool(self.simulation_params.cpus)
        input_pool = [[f, snapshot, self.simulation_params.residue, self.simulation_params.topology] for f in files_out]
        all_coords = pool.map(_extract_coords, input_pool)

        return all_coords

    def rename_folder(self):

        index, name = self.folder.split("_")
        return "{}_Selection".format(index)


@dataclass
class LowestEnergy5(Selection):
    """
    Choose 5% lowest binding energy poses as input for the next simulation.
    Use to select inputs for Rescoring after InducedFitExhaustive and InducedFitFast.
    """
    simulation_params: pv.EnviroBuilder
    folder: str

    def run(self):
        self.folder = self.rename_folder()
        files_out, output_energy = self.extract_poses(percentage=0.05)
        self.choose_refinement_input(files_out, output_energy)
        self.copy_files()
        self.set_next_step()
        return self.simulation_params

    def choose_refinement_input(self, files_out, output_energy):

        n_inputs = self.simulation_params.cpus - 1
        self.directory = os.path.join(os.path.dirname(self.simulation_params.pele_dir), self.folder)

        if len(files_out) > n_inputs:  # pick lowest energy poses to fill all CPU slots
            self.inputs = self.cluster_poses(files_out, output_energy)
        elif len(files_out) < 1:  # if top 5% happen to be less than 1 (e.g. when running tests)
            self.inputs, _ = self.extract_poses(percentage=1, n_poses=1)
        else:  # if number of CPUs matches the number of files in top 5%
            self.inputs = files_out

    def cluster_poses(self, files_out, output_energy):

        # extract coordinates from poses
        all_coords = self.extract_all_coords(files_out)

        # run Gaussian Mixture
        cluster = GaussianMixture(self.simulation_params.cpus-1, covariance_type='full', random_state=42)
        labels = cluster.fit_predict(all_coords)

        # get lowest energy representative of each cluster
        clustered_lig = pd.DataFrame(list(zip(files_out, labels, output_energy)),
                                     columns=["file_name", "cluster_ID", "binding_energy"])
        clustered_lig = clustered_lig.sort_values("binding_energy", ascending=True).groupby("cluster_ID").first()

        output_files = clustered_lig["file_name"].values.tolist()

        return output_files


@dataclass
class Scatter6(Selection):
    """
    Choose input for refinement simulation after the first stage of Allosteric, GPCR and out_in simulations.
    Scan top 75% binding energies, pick n best ones as long as ligand COMs are >= 6 A away from each other.
    """
    simulation_params: pv.EnviroBuilder
    folder: str

    def run(self):
        self.folder = self.rename_folder()
        files_out, output_energy = self.extract_poses(percentage=0.75)
        self.choose_refinement_input(files_out, output_energy)
        self.copy_files()
        self.set_next_step()
        return self.simulation_params

    def choose_refinement_input(self, files_out, output_energy):

        distance = 6
        all_coords = self.extract_all_coords(files_out)
        coords = [list(c[0:3]) for c in all_coords]
        files_out = [os.path.join(self.simulation_params.pele_dir, "results", f) for f in files_out]
        dataframe = pd.DataFrame(list(zip(files_out, output_energy, coords)),
                                 columns=["File", "Binding energy", "1st atom coordinates"])
        dataframe = dataframe.sort_values(["Binding energy"], ascending=True)

        self.inputs = self._check_ligand_distances(dataframe, distance)
        self.directory = os.path.join(os.path.dirname(self.simulation_params.pele_dir),
                                      self.folder)

    def _check_ligand_distances(self, dataframe, distance):

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
                    distances_bool = [d > distance for d in distances]
                    if all(distances_bool):
                        inputs.append(f)
                        input_coords.append(c)
            break  # make sure it stops after running out of files to check

        return inputs
