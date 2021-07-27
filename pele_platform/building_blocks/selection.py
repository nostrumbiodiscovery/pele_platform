from abc import abstractmethod
from dataclasses import dataclass
import glob
import numpy as np
import os
import pandas as pd
import re
from sklearn.mixture import GaussianMixture

from pele_platform.building_blocks import blocks
from pele_platform.Utilities.Helpers import helpers, bestStructs
from pele_platform.Utilities.Parameters.parameters import ParametersBuilder


class Selection(blocks.Block):

    def run(self):
        pass

@dataclass
class Selection(blocks.Block):
    """
    Base class to handle all input selection algorithms, copy files, set next_step, etc.
    """

    simulation_params: ParametersBuilder
    options: dict

    def copy_files(self):
        """
        Copies files selected as self.inputs into a Selection directory.
        """
        self.directory = os.path.join(
            os.path.dirname(self.simulation_params.pele_dir), self.folder_name
        )

        if not os.path.isdir(self.directory):
            os.makedirs(self.directory, exist_ok=True)
        for i in self.inputs:
            os.system("cp {} {}/.".format(i, self.directory))

    def set_next_step(self):
        """
        Sets self.next step so that the inputs can be used by the next Simulation block.
        """
        self.simulation_params.next_step = os.path.join(self.directory, "*.pdb")

    def extract_poses(self, percentage, n_poses=None):
        """
        Extracts n lowest binding energy poses.
        """
        simulation_path = os.path.join(
            self.simulation_params.pele_dir, self.simulation_params.output
        )
        n_best_poses = (
            n_poses
            if n_poses
            else int(
                self.simulation_params.iterations
                * self.simulation_params.pele_steps
                * (self.simulation_params.cpus - 1)
                * percentage
            )
        )

        with helpers.cd(simulation_path):
            files_out, _, _, _, output_energy = bestStructs.main(
                str(self.simulation_params.be_column),
                n_structs=n_best_poses,
                path="",
                topology=self.simulation_params.topology,
                logger=self.simulation_params.logger,
            )

        files_out = [
            os.path.join(self.simulation_params.pele_dir, self.simulation_params.output, f)
            for f in files_out
        ]
        return files_out, output_energy

    def extract_all_coords(self, files_out):
        """
        Extracts ligand coordinates from the previously selected lowest binding energy poses.
        """
        print("Starting coord extraction.")
        snapshot = 0
        input_pool = [
            [
                f,
                snapshot,
                self.simulation_params.residue,
                self.simulation_params.topology,
            ]
            for f in files_out
        ]

        all_coords = helpers.parallelize(_extract_coords, input_pool, 1)
        print("Finished extracting coordinates!")
        return all_coords

    def gaussian_mixture(self, files_out, output_energy):
        # extract coordinates from poses
        all_coords = self.extract_all_coords(files_out)

        # run Gaussian Mixture
        cluster = GaussianMixture(
            self.simulation_params.cpus - 1, covariance_type="full", random_state=42
        )
        labels = cluster.fit_predict(all_coords)

        # get lowest energy representative of each cluster
        clustered_lig = pd.DataFrame(
            list(zip(files_out, labels, output_energy)),
            columns=["file_name", "cluster_ID", "binding_energy"],
        )
        clustered_lig = (
            clustered_lig.sort_values("binding_energy", ascending=True)
            .groupby("cluster_ID")
            .first()
        )

        # save CSV file
        csv_location = os.path.join(
            self.simulation_params.pele_dir, "output", "clustering_output.csv"
        )
        clustered_lig.to_csv(csv_location)

        return clustered_lig["file_name"].values.tolist()

    def rename_folder(self):
        user_folder = self.options.get("working_folder", None) if self.options else None

        if not user_folder:
            index, name = self.folder_name.split("_")
            self.folder_name = "{}_Selection".format(index)
        else:
            self.folder_name = user_folder

        self.simulation_params.folder_name = (
            self.folder_name
        )  # to make it consistent with simulation BBs

    def set_optional_params(self):
        """
        Sets optional params provided by the user in nested yaml (when using 'workflow' flag).
        """
        if self.options:
            for key, value in self.options.items():
                setattr(self.simulation_params, key, value)

    @abstractmethod
    def get_inputs(self):
        """
        Method that will set selected files as self.inputs, ensuring the number of self.inputs is not greater than
        self.n_inputs.
        """
        pass

    def run(self):
        """
        Runs the whole Selection block.
        """
        self.n_inputs = self.simulation_params.cpus - 1
        self.set_optional_params()
        self.rename_folder()
        self.get_inputs()
        self.copy_files()
        self.set_next_step()
        return self.simulation_params


@dataclass
class LowestEnergy5(Selection):
    """
    Choose 5% lowest binding energy poses as input for the next simulation.
    Use to select inputs for Rescoring after LocalExplorationExhaustive and LocalExplorationFast.
    """

    simulation_params: ParametersBuilder
    options: dict
    folder_name: str

    def get_inputs(self):
        files_out, output_energy = self.extract_poses(percentage=0.05)
        self.choose_refinement_input(files_out, output_energy)

    def choose_refinement_input(self, files_out, output_energy):
        if (
            len(files_out) > self.n_inputs
        ):  # pick lowest energy poses to fill all CPU slots
            self.inputs = self.gaussian_mixture(files_out, output_energy)
        elif (
            len(files_out) < 1
        ):  # if top 5% happen to be less than 1 (e.g. when running tests)
            self.inputs, _ = self.extract_poses(percentage=1, n_poses=1)
        else:  # if number of CPUs matches the number of files in top 5%
            self.inputs = files_out


@dataclass
class GMM(Selection):
    """
    Perform Gaussian Mixture (full covariance) clustering on best binding energy poses.
    """

    simulation_params: ParametersBuilder
    options: dict
    folder_name: str

    def get_inputs(self):
        files_out, ou self.gaussian_mixture(files_out, tput_energy = self.extract_poses(percentage=0.05, n_poses=1000)
        self.inputs =output_energy)


@dataclass
class Clusters(Selection):
    """
    Select cluster representatives from 'results' folder as input for the next simulation. If there are more inputs
    than available CPUs, choose the ones with the lowest binding energy.
    """

    simulation_params: ParametersBuilder
    options: dict
    folder_name: str

    def get_inputs(self):
        clusters_dir = os.path.join(
            self.simulation_params.pele_dir, "results/clusters/*.pdb"
        )
        clusters_files = glob.glob(clusters_dir)

        if len(clusters_files) > self.n_inputs:
            self.inputs = self.filter_energies(clusters_files)
        elif len(clusters_files) == 0:
            raise OSError("No files found in {}".format(clusters_dir))
        else:
            self.inputs = clusters_files

    def filter_energies(self, files):
        pattern = r"BindingEnergy([\-0-9]+\.[0-9]+)"
        energies = [float(re.findall(pattern, file)[0]) for file in files]
        zipped = zip(files, energies)
        zipped_list = sorted(zipped, key=lambda x: x[1])
        output = [f for f, e in zipped_list[0 : self.n_inputs]]
        return output


@dataclass
class ScatterN(Selection):
    """
    Choose input for refinement simulation after the first stage of Allosteric, GPCR and out_in simulations.
    Scan top 75% binding energies, pick n best ones as long as ligand COMs are >= 6 A away from each other.
    """

    simulation_params: ParametersBuilder
    options: dict
    folder_name: str
    print("Scatter initialised.")

    def get_inputs(self):
        files_out, output_energy = self.extract_poses(percentage=0.75)
        self.choose_refinement_input(files_out, output_energy)

    def choose_refinement_input(self, files_out, output_energy):
        distance = self.simulation_params.distance
        all_coords = self.extract_all_coords(files_out)
        coords = [list(c[0:3]) for c in all_coords]
        # files_out = [
        #     os.path.join(self.simulation_params.pele_dir, "results", f)
        #     for f in files_out
        # ]
        dataframe = pd.DataFrame(
            list(zip(files_out, output_energy, coords)),
            columns=["File", "Binding energy", "1st atom coordinates"],
        )
        dataframe['Binding energy'] = dataframe['Binding energy'].map(lambda x: float(x))
        dataframe = dataframe.sort_values(["Binding energy"], ascending=True)
        dataframe.to_csv("scatter6_data.csv")
        self.inputs = self._check_ligand_distances(dataframe, distance)

    def _check_ligand_distances(self, dataframe, distance):
        inputs = []
        input_coords = []

        for file, coord in zip(dataframe["File"], dataframe["1st atom coordinates"]):

            if len(inputs) == self.n_inputs:  # get out of the loop, if we have enough inputs already
                break

            if not input_coords:  # first loop
                inputs.append(file)
                input_coords.append(coord)

            else:
                distances = []
                for ic in input_coords:
                    distances.append(
                        abs(np.linalg.norm(np.array(coord) - np.array(ic))))
                print("file: {}, 1st atom coord {}, n of distances {}".format(file, coord, len(distances)))

                distances_bool = [d > distance for d in distances]
                if all(distances_bool):
                    inputs.append(file)
                    input_coords.append(coord)

        return inputs