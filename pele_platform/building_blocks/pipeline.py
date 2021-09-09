from typing import List

from pele_platform.building_blocks.selection import *  # TODO: Is there a better way to do this?
from pele_platform.building_blocks.simulation import *
from pele_platform.context import context


class Pipeline:

    def __init__(self, steps, refinement_steps=None):
        """
        Initializes Pipeline class.

        Parameters
        ----------
        steps : List[dict]
            List of dictionaries containing workflow steps to run if skip_refinement is enabled.
        refinement_steps : List[dict]
            Remaining list of steps
        """
        self.steps = steps
        self.refinement_steps = refinement_steps if refinement_steps else []

        if not context.yaml_parser.skip_refinement:
            self.iterable = self.steps + self.refinement_steps
        else:
            self.iterable = self.steps

    @classmethod
    def make_pipeline(cls, steps, refinement_steps=None):
        """
        Initializes Pipeline class with predefined steps.

        Parameters
        ----------
        steps : List[dict]
            List of dictionaries containing workflow steps to run if skip_refinement is enabled.
        refinement_steps : List[dict]
            Remaining list of steps
        """
        return Pipeline(steps, refinement_steps)

    def run(self):
        """
        Checks the content of the Pipelines and adjusts it, if running in debug mode, then executes everything.

        Returns
        -------
        A list of ParametersBuilder objects containing simulation parameters for each block.
        """
        self.iterable = PipelineChecker(self.iterable).check()
        return self.launch_pipeline()

    def launch_pipeline(self):
        """
        Executes building blocks in the list one by one, appending a deepcopy of Parameters to output after each step.

        Returns
        -------
            A list of Parameters objects
        """
        output = []
        for step in self.iterable:

            simulation_name = step.get("type")
            folder = f"{self.iterable.index(step) + 1}_{simulation_name}"

            simulation = eval(simulation_name)
            options = step.get("options", {})

            simulation(options, folder).run()
            output.append(context.parameters_copy)

        return output


class PipelineChecker:

    def __init__(self, iterable):
        """
        Initializes PipelineChecker class.

        Parameters
        ----------
        iterable : List[dict]
            List of dictionaries containing the simulation type and optional parameters.
        """
        self.iterable = iterable

    def check(self):
        """
        Runs all checks.

        Returns
        -------
            Returns adjusted Pipeline object.
        """

        self.check_pipeline()
        self.check_debug()

        return self.iterable

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
        if context.yaml_parser.debug:
            self.iterable = [self.iterable[0]]  # 1st simulation in the Pipeline only
