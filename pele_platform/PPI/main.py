from pele_platform.Utilities.BuildingBlocks.pipeline import Pipeline
from pele_platform.Utilities.BuildingBlocks.selection import GMM
from pele_platform.Utilities.BuildingBlocks.blocks import Rescoring, InducedFitExhaustive


class PPILauncher:

    def __init__(self, env):
        self.env = env

    def run(self):
        self.env.package = "ppi"

        if not self.env.initial_args.skip_refinement:
            result = Pipeline([InducedFitExhaustive, GMM, Rescoring], self.env).run()
        else:
            result = Pipeline([InducedFitExhaustive], self.env).run()
        return result
