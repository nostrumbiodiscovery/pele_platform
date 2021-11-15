from copy import deepcopy
import glob
import numpy as np
import os
import pandas as pd
import shutil

from pele_platform.Utilities.Helpers.yaml_parser import YamlParser
from pele_platform.Utilities.Parameters.parameters import Parameters


class Equilibrator:

    def __init__(self, args: YamlParser, parameters: Parameters):
        """
        Instantiates the Equilibrator class.

        Parameters
        ----------
        args : YamlParser
            YamlParser object with user-defined args.
        parameters : Parameters
            Parameters object with user-defined args.
        """
        self.args = args
        self.parameters = parameters
        self.equilibration_parameters = None

    def run(self):
        """
        Runs the equilibration workflow to extract cluster values based on protein-ligand contacts and adjust the main
        simulation args accordingly.
        """
        self.generate_equilibration_parameters()
        self.run_equilibration()

        cluster_conditions = self.extract_cluster_conditions(
            os.path.join(
                self.equilibration_parameters.pele_dir,
                self.equilibration_parameters.output,
            )
        )

        return cluster_conditions

    def generate_equilibration_parameters(self):
        """
        Generates args for the short equilibration simulation.
        """
        equilibration_parameters = deepcopy(self.args)
        equilibration_parameters.iterations = 1
        equilibration_parameters.pele_steps = 5
        equilibration_parameters.cluster_conditions = None
        equilibration_parameters.folder = os.path.join(
            self.parameters.pele_dir, self.parameters.output, "preequilibration"
        )
        equilibration_parameters.analyse = False
        equilibration_parameters.no_ppp = True
        equilibration_parameters.randomize = False

        if self.parameters.input:  # in case the original package had multiple inputs or was randomized
            equilibration_parameters.input = self.parameters.input
        else:
            equilibration_parameters.system = self.parameters.system

        self.equilibration_parameters = equilibration_parameters

    def run_equilibration(self):
        """
        Launches the simulation with equilibration args.
        """
        from pele_platform.Adaptive.simulation import run_adaptive

        self.equilibration_parameters = run_adaptive(self.equilibration_parameters)

    def extract_cluster_conditions(self, output_path):
        """
        Extracts contacts and thresholds identified by AdaptivePELE to select the best values for Adaptive clustering.

        Returns
        -------
        cluster_conditions : List[float]
            List of cluster conditions.
        """
        clustering_file = os.path.join(output_path, "0", "clustering", "summary.txt")
        n_cluster_values = len(eval(self.parameters.cluster_values))   # cluster_values are a string

        if not os.path.exists(clustering_file):
            raise FileNotFoundError(
                f"Could not locate {clustering_file} necessary to extract cluster values."
            )

        with open(clustering_file, "r") as file:
            df = pd.read_csv(file, delimiter=r"\s+")
        contacts = df["contacts"].tolist()

        # To avoid cases where there is only one cluster and linspace will fail
        if min(contacts) < max(contacts):
            minimum = min(contacts)
        else:
            minimum = max(max(contacts) - 1, 0)

        # To ensure the final values are not negative
        if minimum < 0:
            minimum = 0

        cluster_conditions = list(
            np.linspace(max(contacts), minimum, num=n_cluster_values - 1)
        )

        print(f"Setting cluster conditions {cluster_conditions}.")
        return cluster_conditions

    def set_next_step(self):
        """
        Extracts the best poses generated during the equilibration, copies them over to the main simulation inputs dir
        and overwrites the system argument.
        """

        from pele_platform.analysis.analysis import Analysis

        next_inputs_path = os.path.join(
            self.equilibration_parameters.pele_dir, "selected_poses"
        )
        analysis = Analysis.from_parameters(self.equilibration_parameters)
        analysis.generate_top_poses(
            next_inputs_path, self.equilibration_parameters.cpus - 1
        )

        inputs_dir = os.path.join(self.parameters.pele_dir, "input")
        old_systems = glob.glob(
            os.path.join(
                inputs_dir,
                os.path.basename(self.args.system.replace(".pdb", "*")),
            )
        )

        for old_system in old_systems:
            os.remove(old_system)

        inputs = glob.glob(os.path.join(next_inputs_path, "*.pdb"))

        for input_file in inputs:
            shutil.copy(input_file, inputs_dir)

        return [os.path.basename(input_file) for input_file in inputs]
