from copy import deepcopy
from dataclasses import dataclass

from pele_platform.Utilities.Parameters import pele_env as pv


@dataclass
class Pipeline:
    iterable: list
    env: pv.EnviroBuilder

    def run(self):
        output = []
        for simulation in self.iterable:
            folder = "{}_{}".format(self.iterable.index(simulation) + 1, simulation.__name__)
            self.env = simulation(self.env, folder).run()
            output.append(deepcopy(self.env))

        return output
