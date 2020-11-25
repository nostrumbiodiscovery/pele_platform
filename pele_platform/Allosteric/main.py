from dataclasses import dataclass
import os

from pele_platform.Utilities.BuildingBlocks import blocks as bb
import pele_platform.Utilities.Parameters.pele_env as pv


class AllostericLauncher:

    def __init__(self, env):
        self.env = env

    def run(self):
        """
        Launch allosteric simulation.
        1) Run global exploration to identify the most important pockets
        2) Run induced fit simulation to find deep pockets
        """
        self.env.package = "allosteric"
        result = Pipeline([bb.GlobalExploration, bb.InducedFitExhaustive], self.env)

        return result

def Pipeline(iterable, env):
    for simulation in iterable:
        env = simulation(env).run()
