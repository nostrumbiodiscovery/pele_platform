from abc import abstractmethod
from dataclasses import dataclass
import glob
import numpy as np
import os
import pandas as pd
import shutil

from pele_platform.building_blocks import blocks
import pele_platform.Utilities.Parameters.parameters as pv
from pele_platform.analysis.analysis import Analysis


class Selection(blocks.Block):
    """
    Base class to handle all input selection algorithms, copy files, set next_step, etc.
    """

    def __init__(
        self,
        parameters_builder: pv.ParametersBuilder,
        options: dict,
        folder_name: str,
        env: pv.Parameters,
    ):
        self.options = options
        self.folder_name = folder_name
        self.builder = parameters_builder
        self.env = env
        self.n_inputs = self.env.cpus - 1
        self.inputs = None
        self.analysis = Analysis.from_parameters(self.env)

    def copy_files(self):
        """
        Copies files selected as self.inputs into a Selection directory.
        """
        if not os.path.isdir(self.env.pele_dir):
            os.makedirs(self.env.pele_dir, exist_ok=True)

        for i in self.inputs:
            try:
                shutil.copy(i, self.env.pele_dir)
            except shutil.SameFileError:
                pass

    def set_next_step(self):
        """
        Sets self.next step so that the inputs can be used by the next Simulation block.
        """
        self.env.next_step = os.path.join(self.env.pele_dir, "*.pdb")

    def rename_folder(self):
        """
        Check if user specified a custom folder name for this simulation block. If not, use the automatically
        generated one.
        """
        if self.options:
            user_folder = self.options.get("working_folder", None)
            self.env.folder_name = user_folder if user_folder else self.folder_name
        else:
            self.env.folder_name = self.folder_name

        if self.env.pele_dir == self.builder.pele_dir:
            self.env.pele_dir = os.path.join(self.builder.pele_dir, self.env.folder_name)
        else:
            self.env.pele_dir = os.path.join(os.path.dirname(self.env.pele_dir), self.env.folder_name)

        self.env.inputs_dir = os.path.join(self.env.pele_dir, "input")

    # def rename_folder(self):
    #     user_folder = self.options.get("working_folder", None) if self.options else None
    #
    #     if not user_folder:
    #         index, name = self.folder_name.split("_")
    #         self.folder_name = "{}_Selection".format(index)
    #     else:
    #         self.folder_name = user_folder
    #
    #     self.env.folder_name = (
    #         self.folder_name
    #     )  # to make it consistent with simulation BBs

    def set_optional_params(self):
        """
        Sets optional params provided by the user in nested yaml (when using 'workflow' flag).
        """
        if self.options:
            for key, value in self.options.items():
                setattr(self.env, key, value)

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
        self.set_optional_params()
        self.get_inputs()
        self.rename_folder()
        self.copy_files()
        self.set_next_step()
        return self.builder, self.env


@dataclass
class LowestEnergy(Selection):
    """
    Choose lowest binding energy poses as input for the next simulation.
    Use to select inputs for Rescoring after LocalExplorationExhaustive and LocalExplorationFast.
    """

    def __init__(
        self,
        parameters_builder: pv.ParametersBuilder,
        options: dict,
        folder_name: str,
        env: pv.Parameters,
    ):
        super().__init__(parameters_builder, options, folder_name, env)

    def get_inputs(self):
        top_poses_folder = os.path.join(self.env.pele_dir, "temp")

        self.analysis.generate_top_poses(top_poses_folder, n_poses=self.n_inputs)
        self.inputs = glob.glob(os.path.join(top_poses_folder, "*.pdb"))


@dataclass
class GMM(Selection):
    """
    Perform Gaussian Mixture (full covariance) clustering on best binding energy poses.
    """

    def __init__(
        self,
        parameters_builder: pv.ParametersBuilder,
        options: dict,
        folder_name: str,
        env: pv.Parameters,
    ):
        super().__init__(parameters_builder, options, folder_name, env)

    def get_inputs(self):
        temp_dir = "temp"

        self.analysis.generate_clusters(
            temp_dir,
            analysis_nclust=self.n_inputs,
            clustering_type="GaussianMixture",
        )
        self.inputs = glob.glob(
            os.path.join(temp_dir, "cluster*.pdb")
        )


@dataclass
class Clusters(Selection):
    """
    Select cluster representatives from 'results' folder as input for the next simulation. If there are more inputs
    than available CPUs, choose the ones with the lowest binding energy.
    """

    def __init__(
        self,
        parameters_builder: pv.ParametersBuilder,
        options: dict,
        folder_name: str,
        env: pv.Parameters,
    ):
        super().__init__(parameters_builder, options, folder_name, env)
        self.inputs = None

    def get_inputs(self):

        clusters_dir = os.path.join(self.env.pele_dir, "results/clusters/cluster*.pdb")
        clusters_files = glob.glob(clusters_dir)

        if len(clusters_files) > self.n_inputs:
            df = pd.read_csv(
                os.path.join(os.path.dirname(clusters_dir), "top_selections.csv")
            )
            df = df.nsmallest(n=self.n_inputs, columns="Binding Energy")
            files_to_select = [
                os.path.join(os.path.dirname(clusters_dir), f"cluster_{label}.pdb")
                for label in df["Cluster label"]
            ]

            self.inputs = files_to_select
        else:
            self.inputs = clusters_files


@dataclass
class ScatterN(Selection):
    """
    Choose input for refinement simulation after the first stage of Allosteric, GPCR and out_in simulations.
    Scan top 75% binding energies, pick n best ones as long as ligand COMs are >= 6 A away from each other.
    """

    def __init__(
        self,
        parameters_builder: pv.ParametersBuilder,
        options: dict,
        folder_name: str,
        env: pv.Parameters,
    ):
        super().__init__(parameters_builder, options, folder_name, env)
        self.inputs = None

    def get_inputs(self):
        pass

    def _check_ligand_distances(self, dataframe, distance):
        inputs = []
        input_coords = []

        for file, coord in zip(dataframe["File"], dataframe["1st atom coordinates"]):

            if (
                len(inputs) == self.n_inputs
            ):  # get out of the loop, if we have enough inputs already
                break

            if not input_coords:  # first loop
                inputs.append(file)
                input_coords.append(coord)

            else:
                distances = []
                for ic in input_coords:
                    distances.append(
                        abs(np.linalg.norm(np.array(coord) - np.array(ic)))
                    )
                print(
                    "file: {}, 1st atom coord {}, n of distances {}".format(
                        file, coord, len(distances)
                    )
                )

                distances_bool = [d > distance for d in distances]
                if all(distances_bool):
                    inputs.append(file)
                    input_coords.append(coord)

        return inputs
