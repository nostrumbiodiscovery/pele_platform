from copy import deepcopy
from dataclasses import dataclass
from pele_platform.Utilities.Parameters import pele_env as pv
import pele_platform.Errors.custom_errors as ce


@dataclass
class Pipeline:
    iterable: list
    env: pv.EnviroBuilder

    def run(self):
        self.check_pipeline()
        self.check_debug()
        output = []
        for simulation in self.iterable:
            folder = "{}_{}".format(self.iterable.index(simulation) + 1, simulation.__name__)
            self.env = simulation(self.env, folder).run()
            output.append(deepcopy(self.env))
        return output


    def check_pipeline(self):
        
        block_types = [block.__bases__[-1].__name__ for block in self.iterable]
        
        if len(block_types) == 0:
            raise ce.PipelineError("Pipeline doesn't contain any BuildingBlocks.")

        if not block_types[0] == "Simulation" or not block_types[-1] == "Simulation":
            raise ce.PipelineError("Pipeline should begin and end with a Simulation block, e.g. GlobalExploration.")

        for index, block in enumerate(block_types):
            if (index % 2 == 0 and block != "Simulation") or (index % 2 == 1 and block != "Selection"):
                raise ce.PipelineError("There is something wrong with your Pipeline.\n1. It should begin and end with a Simulation block, e.g. GlobalExploration.\n2. Simulation blocks should be separated by Selection blocks.")
            

    def check_debug(self):
        """
        If input.yaml contains debug: true flag, the execution of the Pipeline should stop after setting up
        parameters and folders for the first simulation.
        """
        if self.env.initial_args.debug:
            self.env.debug = self.env.initial_args.debug
            self.iterable = [self.iterable[0]]  # 1st simulation in the Pipeline only
