from pele_platform.Utilities.BuildingBlocks.pipeline import Pipeline
from pele_platform.Utilities.BuildingBlocks.selection import LowestEnergy5
from pele_platform.Utilities.BuildingBlocks.blocks import InducedFitExhaustive, InducedFitFast, Rescoring


class InducedFitExhaustiveLauncher:

    def __init__(self, env):
        self.env = env

    def run(self):
        self.env.package = "induced_fit_exhaustive"

        result = Pipeline([InducedFitExhaustive, LowestEnergy5, Rescoring], self.env).run()
        return result


class InducedFitFastLauncher:

    def __init__(self, env):
        self.env = env

    def run(self):
        self.env.package = "induced_fit_fast"

        result = Pipeline([InducedFitFast, LowestEnergy5, Rescoring], self.env).run()
        return result
