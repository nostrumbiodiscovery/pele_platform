from pele_platform.Utilities.BuildingBlocks.pipeline import Pipeline
from pele_platform.Utilities.BuildingBlocks.selection import Scatter6
from pele_platform.Utilities.BuildingBlocks.blocks import GPCR, Rescoring


class GPCRLauncher:

    def __init__(self, env):
        self.env = env

    def run(self):
        self.env.package = "gpcr_orth"

        result = Pipeline([GPCR, Scatter6, Rescoring], self.env).run()
        return result
