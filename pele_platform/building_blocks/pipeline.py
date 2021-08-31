from copy import deepcopy
from typing import List

from pele_platform.building_blocks.selection import *
from pele_platform.building_blocks.simulation import *
from pele_platform.Utilities.Parameters.parameters import ParametersBuilder


class Pipeline:
    def __init__(self, iterable: List[dict], builder: ParametersBuilder):
        """
        Initializes Pipeline class.

        Parameters
        ----------
        iterable : List[dict]
            List of dictionaries containing workflow steps and their optional parameters.
        builder : ParametersBuilder
            ParametersBuilder object.
        """
        self.iterable = iterable
        self.builder = builder
        self.env = None

    def run(self):
        """
        Run the whole Pipeline.

        Returns
        -------
        A list of ParametersBuilder objects containing simulation parameters for each block.
        """
        self.check_pipeline()
        self.check_debug()
        return self.launch_pipeline()
        # TODO: Separate out the checker

    def launch_pipeline(self):
        output = []
        for step in self.iterable:
            simulation_name = step.get("type")
            simulation = eval(simulation_name)
            options = step.get("options", {})
            folder = "{}_{}".format(self.iterable.index(step) + 1, simulation_name)
            print("AAAAAAAAAAAAAAAA folder in Pipeline", folder)
            self.builder, self.env = simulation(self.builder, options, folder, self.env).run()
            output.append(deepcopy(self.env))
        return output

    def check_pipeline(self):
        """
        Checks the pipeline to ensure it's not empty and the simulation are selection blocks occur alternately.
        """

        if len(self.iterable) == 0:
            raise custom_errors.PipelineError("Pipeline doesn't contain any building_blocks.")

        for element in self.iterable:
            if "type" not in element.keys():
                raise custom_errors.PipelineError(
                    "It seems that you forgot to specify simulation type in one of the workflow "
                    "steps, for example - type: 'LocalExplorationExhaustive'"
                )
            if element["type"] == "ScatterN":
                dist = element.get("options", {}).get("distance", None)
                if not (isinstance(dist, float) or isinstance(dist, int)):
                    raise ValueError(
                        "ScatterN requires optional parameter distance which can be either an integer or "
                        "a float. Please refer to the PELE Platform documentation for more details on "
                        "building_blocks."
                    )

        class_names = [eval(element["type"]) for element in self.iterable]
        block_types = [block.__bases__[-1].__name__ for block in class_names]

        if not block_types[0] == "Simulation" or not block_types[-1] == "Simulation":
            raise custom_errors.PipelineError(
                "Workflow should begin and end with a Simulation block, e.g. GlobalExploration."
            )

        for index, block in enumerate(block_types):
            if (index % 2 == 0 and block != "Simulation") or (
                index % 2 == 1 and block != "Selection"
            ):
                raise custom_errors.PipelineError(
                    "There is something wrong with your Workflow.\n1. It should begin and end with a Simulation "
                    "block, e.g. GlobalExploration.\n2. Simulation blocks should be separated by Selection blocks."
                )

    def check_debug(self):
        """
        If input.yaml contains debug: true flag, the execution of the Pipeline should stop after setting up
        parameters and folders for the first simulation.
        """
        if self.builder.initial_args.debug:
            self.iterable = [self.iterable[0]]  # 1st simulation in the Pipeline only
