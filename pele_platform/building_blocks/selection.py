from abc import abstractmethod
from dataclasses import dataclass
import glob
import numpy as np
import os
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
        self.simulation_params = env
        self.directory = os.path.join(os.path.dirname(self.simulation_params.pele_dir), self.folder_name)
        print("AAAAAAAAAAAAAAA self.directory in Selection init", self.directory)
        self.n_inputs = self.simulation_params.cpus - 1
        self.inputs = None
        self.analysis = Analysis.from_parameters(self.simulation_params)

    def copy_files(self):
        """
        Copies files selected as self.inputs into a Selection directory.
        """
        if not os.path.isdir(self.directory):
            os.makedirs(self.directory, exist_ok=True)

        for i in self.inputs:
            try:
                shutil.copy(i, self.directory)
            except shutil.SameFileError:
                pass

    def set_next_step(self):
        """
        Sets self.next step so that the inputs can be used by the next Simulation block.
        """
        self.simulation_params.next_step = os.path.join(self.directory, "*.pdb")

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
        self.set_optional_params()
        self.rename_folder()
        self.get_inputs()
        self.copy_files()
        self.set_next_step()
        return self.builder, self.simulation_params


@dataclass
class LowestEnergy(Selection):
    """
    Choose lowest binding energy poses as input for the next simulation.
    Use to select inputs for Rescoring after LocalExplorationExhaustive and LocalExplorationFast.
    """

    def __init__(self, parameters_builder: pv.ParametersBuilder, options: dict, folder_name: str, env: pv.Parameters):
        super().__init__(parameters_builder, options, folder_name, env)

    def get_inputs(self):
        self.analysis.generate_top_poses(os.path.join(self.simulation_params.pele_dir, self.folder_name),
                                         n_poses=self.n_inputs)
        self.inputs = glob.glob(os.path.join(self.simulation_params.folder_name, "cluster*.pdb"))


@dataclass
class GMM(Selection):
    """
    Perform Gaussian Mixture (full covariance) clustering on best binding energy poses.
    """

    def __init__(self, parameters_builder: pv.ParametersBuilder, options: dict, folder_name: str, env: pv.Parameters):
        super().__init__(parameters_builder, options, folder_name, env)

    def get_inputs(self):
        self.analysis.generate_clusters(self.simulation_params.folder_name, analysis_nclust=self.n_inputs,
                                        clustering_type="gmm")


@dataclass
class Clusters(Selection):
    """
    Select cluster representatives from 'results' folder as input for the next simulation. If there are more inputs
    than available CPUs, choose the ones with the lowest binding energy.
    """

    def __init__(self, parameters_builder: pv.ParametersBuilder, options: dict, folder_name: str, env: pv.Parameters):
        super().__init__(parameters_builder, options, folder_name, env)
        self.inputs = None

    def get_inputs(self):
        clusters_dir = os.path.join(self.simulation_params.pele_dir, "results/clusters/cluster*.pdb")
        clusters_files = glob.glob(clusters_dir)
        self.inputs = clusters_files


@dataclass
class ScatterN(Selection):
    """
    Choose input for refinement simulation after the first stage of Allosteric, GPCR and out_in simulations.
    Scan top 75% binding energies, pick n best ones as long as ligand COMs are >= 6 A away from each other.
    """

    def __init__(self, parameters_builder: pv.ParametersBuilder, options: dict, folder_name: str, env: pv.Parameters):
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
