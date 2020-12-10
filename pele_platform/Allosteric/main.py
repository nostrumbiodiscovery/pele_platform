from pele_platform.Utilities.BuildingBlocks.pipeline import Pipeline
from pele_platform.Utilities.BuildingBlocks.selection import Scatter6
from pele_platform.Utilities.BuildingBlocks.blocks import GlobalExploration, InducedFitExhaustive


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
        if not self.env.initial_args.skip_refinement:
            result = Pipeline([GlobalExploration, Scatter6, InducedFitExhaustive], self.env).run()
        else:
            result = Pipeline([GlobalExploration], self.env).run()

        return result
