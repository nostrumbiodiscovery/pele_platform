from copy import deepcopy
import glob
import numpy as np
import os
import pandas as pd
import shutil

from pele_platform.Utilities.Helpers.yaml_parser import YamlParser


class Equilibrator:

    def __init__(self, parameters: YamlParser, pele_dir, preprocessed_system):
        """
        Instantiates the Equilibrator class.

        Parameters
        ----------
        parameters : YamlParser
            YamlParser object with user-defined parameters.
        pele_dir : str
            Path to the PELE directory in which the main simulation is performed.
        preprocessed_system : str
            Path to the system preprocessed by PPP.
        """
        self.parameters = parameters
        self.pele_dir = pele_dir
        self.equilibration_parameters = None
        self.preprocessed_system = preprocessed_system

    def run(self):
        """
        Runs the equilibration workflow to extract cluster values based on protein-ligand contacts and adjust the main
        simulation parameters accordingly.
        """
        self.generate_equilibration_parameters()
        self.run_equilibration()

        cluster_values, cluster_conditions = self.extract_cluster_values(
            os.path.join(self.equilibration_parameters.pele_dir, self.equilibration_parameters.output))

        return cluster_values, cluster_conditions

    def generate_equilibration_parameters(self):
        """
        Generates parameters for the short equilibration simulation.
        """
        equilibration_parameters = deepcopy(self.parameters)
        equilibration_parameters.iterations = 1
        equilibration_parameters.pele_steps = 5
        equilibration_parameters.cluster_values = None
        equilibration_parameters.folder = os.path.join(self.pele_dir, "QuickEquilibration")
        equilibration_parameters.analyse = False
        equilibration_parameters.no_ppp = True
        equilibration_parameters.system = self.preprocessed_system

        self.equilibration_parameters = equilibration_parameters

    def run_equilibration(self):
        """
        Launches the simulation with equilibration parameters.
        """
        from pele_platform.Adaptive.simulation import run_adaptive

        self.equilibration_parameters = run_adaptive(self.equilibration_parameters)

    @staticmethod
    def extract_cluster_values(output_path):
        """
        Extracts contacts and thresholds identified by AdaptivePELE to select the best values for Adaptive clustering.

        Returns
        -------
        cluster_values : List[float]
            List of cluster values (thresholds) to be used in the actual, user-defined simulation.
        cluster_conditions : List[float]
            List of cluster conditions.
        """
        clustering_file = os.path.join(output_path, "0", "clustering", "summary.txt")

        if not os.path.exists(clustering_file):
            raise FileNotFoundError(f"Could not locate {clustering_file} necessary to extract cluster values.")

        with open(clustering_file, "r") as file:
            df = pd.read_csv(file, delimiter=r"\s+")

        contacts = df["contacts"].tolist()
        threshold = df["threshold"].tolist()

        cluster_values = [min(threshold), np.mean(threshold), max(threshold)]
        cluster_conditions = [max(contacts), np.mean(contacts), min(contacts)]

        print(f"Setting cluster values {cluster_values} and cluster conditions {cluster_conditions}.")

        return cluster_values, cluster_conditions

    def set_next_step(self):
        """
        Extracts the best poses generated during the equilibration, copies them over to the main simulation inputs dir
        and overwrites the system argument.
        """

        from pele_platform.analysis.analysis import Analysis

        next_inputs_path = os.path.join(self.equilibration_parameters.pele_dir, "selected_poses")
        analysis = Analysis.from_parameters(self.equilibration_parameters)
        analysis.generate_top_poses(next_inputs_path, self.equilibration_parameters.cpus - 1)

        inputs_dir = os.path.join(self.pele_dir, "input")
        old_systems = glob.glob(os.path.join(inputs_dir, os.path.basename(self.parameters.system.replace(".pdb", "*"))))

        for old_system in old_systems:
            os.remove(old_system)

        inputs = glob.glob(os.path.join(next_inputs_path, "*.pdb"))

        for input_file in inputs:
            shutil.copy(input_file, inputs_dir)

        return [os.path.basename(input_file) for input_file in inputs]
