from collections import defaultdict
import glob
import numpy as np
import os
import re

from pele_platform.analysis.data import DataHandler
from pele_platform.Utilities.Helpers import helpers
from pele_platform.analysis.clustering import MeanShiftClustering
from pele_platform.analysis.analysis import Analysis


class WaterAnalysis:
    def __init__(
        self,
        simulation_output,
        be_column=4,
        trajectory="trajectory.pdb",
        report="report",
        topology=None,
        skip_initial_structures=True,
        water_ids=None,
    ):
        """
        Initializes WaterAnalysis class.

        Parameters
        -----------
        simulation_output : str
            Path to simulation output folder.
        be_column : int
            Column with energy metric, default 4.
        trajectory : str
            Trajectory name defaults to "trajectory.pdb", but you should use "trajectory.xtc" if using XTC format.
        report : str
            Report name.
        topology : str
            Path to the topology file (mandatory, if using XTC trajectories).
        skip_initial_structures : bool
            Toggle to skip extraction of the very first model.
        water_ids : list[tuple[str, int]]
            The list of water ids to track. Each water id is defined with
            a tuple that contains the PDB chain and the residue number
            corresponding to each water molecule to track.
        """
        self.output = simulation_output
        self.be_column = be_column
        self.trajectory = trajectory
        self.report = report
        self.topology = topology
        self.skip_initial_structures = skip_initial_structures
        self.water_ids = water_ids

        self.data_handler = DataHandler(
            sim_path=self.output,
            report_name=self.report,
            trajectory_name=self.trajectory,
            be_column=self.be_column,
            skip_initial_structures=self.skip_initial_structures,
        )

    @classmethod
    def from_parameters(cls, simulation_parameters):
        """
        Initializes WaterAnalysis class from simulation parameters.

        Parameters
        ----------
        simulation_parameters : Parameters
            Simulation parameters object.

        Returns
        --------
        WaterAnalysis object.
        """
        analysis = WaterAnalysis(
            simulation_output=os.path.join(
                simulation_parameters.pele_dir, simulation_parameters.output
            ),
            be_column=simulation_parameters.be_column,
            trajectory=simulation_parameters.traj_name,
            report=simulation_parameters.report_name,
            topology=simulation_parameters.topology,
            skip_initial_structures=not simulation_parameters.test,
            water_ids=simulation_parameters.water_ids_to_track,
        )

        return analysis

    def analyse_waters(self, path):
        """
        Runs the whole water analysis workflow:
            - Creates folder.
            - Extracts water coordinates (only tracked IDs).
            - Clusters using mean shift algorithm.
            - TODO: Analyses energy and entropy of each water site.
        """
        helpers.check_make_folder(path)

        coordinates = self.extract_water_coordinates(
            skip_initial_structures=self.skip_initial_structures
        )
        clusters, water_sites = self.get_water_sites(coordinates, path)
        print(f"Generated {len(clusters)} clusters.")

    def extract_water_coordinates(self, skip_initial_structures=True):
        """
        Extracts water coordinates from all trajectories.

        Returns
        -------
        """
        trajectories = glob.glob(
            os.path.join(self.output, "*", self.trajectory.replace(".", "*"))
        )

        if self.topology:
            coordinates = self.extract_waters_from_XTC(
                trajectories,
                water_ids=self.water_ids,
                skip_initial_structures=skip_initial_structures,
            )
        else:
            coordinates = self.extract_waters_from_PDB(
                trajectories,
                water_ids=self.water_ids,
                skip_initial_structures=skip_initial_structures,
            )
        print("Water coordinates extracted.")
        return coordinates

    @staticmethod
    def extract_waters_from_PDB(all_files, water_ids, skip_initial_structures=True):
        """
        Extracts coordinates of water molecules from PDB trajectories.

        Parameters
        ----------
        all_files : List[str]
            List of all trajectory files.
        skip_initial_structures : bool
            Toggle to exclude the first model in each trajectory.
        water_ids : list[tuple[str, int]]
            The list of water ids to track. Each water id is defined with
            a tuple that contains the PDB chain and the residue number
            corresponding to each water molecule to track.

        Returns
        -------
            An array of [M, N, D] shape, where M - number of models, N - number of atoms in the molecule (3), D - number
            of dimensions (3).
        """
        water_coordinates = defaultdict(list)

        if len(all_files) == 0:
            raise FileNotFoundError("No output files detected.")

        for f in all_files:
            with open(f, "r") as file:
                file_content = file.read()
                models = re.findall(r"(MODEL.+?ENDMDL)", file_content, re.DOTALL)

                for model in models:
                    model_number = re.search(r"MODEL\s+(\d+)", model).group(1)
                    if model_number == 1 and skip_initial_structures:
                        continue

                    water_lines = re.findall(r".+HOH.+", model)

                    for line in water_lines:
                        chain = line[21]
                        residue_number = line[22:26].strip()

                        if (chain, int(residue_number)) in water_ids:
                            try:
                                x, y, z = (
                                    float(line[30:38]),
                                    float(line[38:46]),
                                    float(line[46:54]),
                                )
                            except ValueError:
                                print(f"Warning: invalid PDB format in trajectory {f}.")

                        point = np.array((x, y, z))
                        water_coordinates[
                            f"{f}_{model_number}_{chain}_{residue_number}"
                        ].append(point)

        all_coordinates = [
            np.array(coordinates) for coordinates in water_coordinates.values()
        ]
        all_coordinates = np.array(all_coordinates)

        n_models, n_atoms, dimensions = all_coordinates.shape
        assert dimensions == 3, "Extracted water coordinates are not 3-dimensional."
        assert (
            n_atoms == 3
        ), "Water molecules do not seem to have the correct number of atoms."

        return all_coordinates

    @staticmethod
    def extract_waters_from_XTC(all_files, water_ids, skip_initial_structures=True):
        """
        Extracts coordinates of water molecules from XTC trajectories.

        Parameters
        ----------
        all_files : List[str]
            List of all trajectory files.
        skip_initial_structures : bool
            Toggle to exclude the first model in each trajectory.
        water_ids : list[tuple[str, int]]
            The list of water ids to track. Each water id is defined with
            a tuple that contains the PDB chain and the residue number
            corresponding to each water molecule to track.

        Returns
        -------
            An array of [M, N, D] shape, where M - number of models, N - number of atoms in the molecule (3), D - number
            of dimensions (3).
        """
        water_coordinates = list()

        # TODO: Implement extraction from XTC trajectories, preferably as a static method.

        return water_coordinates

    @staticmethod
    def get_water_sites(coordinates, path):
        """
        Cluster water coordinates using mean shift algorithm and write the centroids to PDB file.

        Parameters
        ----------
        coordinates : List[array]
            List of numpy arrays with water coordinates.
        path : str
            Path to save the PDB files.
        """
        watersites_data = list()

        clustering = MeanShiftClustering(1.5)  # hard=coded bandwidth value for water
        water_clusters, estimator = clustering.get_clusters(coordinates)

        populations = Analysis._get_cluster_populations(water_clusters)
        output_path = os.path.join(path, "water_sites.pdb")
        Analysis._write_centroids(populations, estimator, output_path)
        centroids = estimator.cluster_centers_

        for cluster, centroid in enumerate(centroids):
            watersites_data.append([cluster, *centroid, populations[cluster]])

        return water_clusters, watersites_data
