from dataclasses import dataclass
import os

from pele_platform.Utilities.BuildingBlocks import blocks as bb
import pele_platform.Utilities.Parameters.pele_env as pv


@dataclass
class AllostericLauncher:
    args: pv.EnviroBuilder

    def run(self):
        """
        Launch allosteric simulation.
        1) Run global exploration to identify the most important pockets
        2) Run induced fit simulation to find deep pockets
        """
        self.global_simulation = bb.GlobalExploration(self.args).run()

        if not self.args.skip_refinement:
            self.refinement_simulation = bb.InducedFitExhaustive(self.args).run()
        else:
            self.refinement_simulation = None

        return self.global_simulation, self.refinement_simulation

