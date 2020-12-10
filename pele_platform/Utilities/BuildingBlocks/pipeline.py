from copy import deepcopy
from dataclasses import dataclass
from pele_platform.Utilities.Parameters import pele_env as pv


@dataclass
class Pipeline:
    iterable: list
    env: pv.EnviroBuilder

    def run(self):
        self.check_debug()
        output = []
        for simulation in self.iterable:
            folder = "{}_{}".format(self.iterable.index(simulation) + 1, simulation.__name__)
            self.env = simulation(self.env, folder).run()
            output.append(deepcopy(self.env))
        return output

    def check_debug(self):
        """
        If input.yaml contains debug: true flag, the execution of the Pipeline should stop after setting up
        parameters and folders for the first simulation.
        """
        if self.env.initial_args.debug:
            self.env.debug = self.env.initial_args.debug
            self.iterable = [self.iterable[0]]  # 1st simulation in the Pipeline only
