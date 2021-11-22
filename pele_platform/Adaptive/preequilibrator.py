from copy import deepcopy
from itertools import chain
import numpy as np
import os
import pandas as pd

from pele_platform.Utilities.Helpers.yaml_parser import YamlParser
from pele_platform.Utilities.Parameters.parameters import Parameters


class PreEquilibrator:

    def __init__(self, args: YamlParser, parameters: Parameters):
        """
        Instantiates the PreEquilibrator class.

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
        equilibration_parameters.box_radius = self.parameters.box_radius
        equilibration_parameters.box_center = self.parameters.box_center
        equilibration_parameters.analyse = False

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
        sizes = df["size"].tolist()

        contacts_distribution = list(chain.from_iterable([[contact] * size for contact, size in zip(contacts, sizes)]))

        maximum = np.quantile(contacts_distribution, 0.9)
        minimum = np.quantile(contacts_distribution, 0.3)

        if (maximum - minimum < 0.8) and (maximum > 1.0):
            minimum = maximum - 0.9

        elif (maximum - minimum < 0.3) and (maximum > 0.5):
            minimum = maximum - 0.3
            
        cluster_conditions = list(
            np.linspace(maximum, minimum, num=n_cluster_values - 1)
        )

        print(f"Setting cluster conditions {cluster_conditions}.")
        return cluster_conditions
