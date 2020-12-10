from pele_platform.Utilities.BuildingBlocks.blocks import OutIn, Rescoring
from pele_platform.Utilities.BuildingBlocks.pipeline import Pipeline
from pele_platform.Utilities.BuildingBlocks.selection import Scatter6


class OutInLauncher:

    def __init__(self, env):
        self.env = env

    def run(self):
        self.env.package = "out_in"

        result = Pipeline([OutIn, Scatter6, Rescoring], self.env).run()
        return result
