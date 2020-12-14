from pele_platform.Utilities.BuildingBlocks.pipeline import Pipeline
from pele_platform.Utilities.BuildingBlocks.selection import Scatter6
from pele_platform.Utilities.BuildingBlocks.preparation import Preparation
from pele_platform.Utilities.BuildingBlocks.blocks import Rescoring, InducedFitExhaustive


class PPILauncher:

    def __init__(self, env):
        self.env = env

    def run(self):
        self.env.package = "ppi"

        if not self.env.initial_args.skip_refinement:
            result = Pipeline([Preparation, InducedFitExhaustive, Scatter6, Rescoring], self.env).run()
        else:
            result = Pipeline([Preparation, InducedFitExhaustive], self.env).run()
        return result
