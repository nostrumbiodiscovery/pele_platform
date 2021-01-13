from copy import deepcopy
from dataclasses import dataclass
import pele_platform.Errors.custom_errors as ce
from pele_platform.Utilities.BuildingBlocks.simulation import *
from pele_platform.Utilities.BuildingBlocks.selection import *
from pele_platform.Utilities.Parameters import pele_env as pv


@dataclass
class Pipeline:
    iterable: list
    env: pv.EnviroBuilder

    def run(self):
        self.check_pipeline()
        self.check_debug()
        output = self.launch_pipeline()
        return output

    def launch_pipeline(self):
        output = []
        for step in self.iterable:
            simulation_name = step.get('type')
            simulation = eval(simulation_name)  # ugly
            options = step.get('options', {})
            folder = "{}_{}".format(self.iterable.index(step) + 1, simulation_name)
            self.env = simulation(self.env, options, folder).run()
            output.append(deepcopy(self.env))
        return output

    def check_pipeline(self):

        if len(self.iterable) == 0:
            raise ce.PipelineError("Pipeline doesn't contain any BuildingBlocks.")

        for element in self.iterable:
            if 'type' not in element.keys():
                raise ce.PipelineError("It seems that you forgot to specify simulation type in one of the workflow "
                                       "steps, for example - type: 'InducedFitExhaustive'")
            if element['type'] == "ScatterN":
                dist = element.get('options', {}).get('distance', None)
                if not(isinstance(dist, float) or isinstance(dist, int)):
                    raise ValueError("ScatterN requires optional parameter distance which can be either an integer or "
                                     "a float. Please refer to the PELE Platform documentation for more details on "
                                     "BuildingBlocks.")

        class_names = [eval(element['type']) for element in self.iterable]
        block_types = [block.__bases__[-1].__name__ for block in class_names]

        if not block_types[0] == "Simulation" or not block_types[-1] == "Simulation":
            raise ce.PipelineError("Workflow should begin and end with a Simulation block, e.g. GlobalExploration.")

        for index, block in enumerate(block_types):
            if (index % 2 == 0 and block != "Simulation") or (index % 2 == 1 and block != "Selection"):
                raise ce.PipelineError(
                    "There is something wrong with your Workflow.\n1. It should begin and end with a Simulation "
                    "block, e.g. GlobalExploration.\n2. Simulation blocks should be separated by Selection blocks.")

    def check_debug(self):
        """
        If input.yaml contains debug: true flag, the execution of the Pipeline should stop after setting up
        parameters and folders for the first simulation.
        """
        if self.env.initial_args.debug:
            self.env.debug = self.env.initial_args.debug
            self.iterable = [self.iterable[0]]  # 1st simulation in the Pipeline only
