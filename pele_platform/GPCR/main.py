import pele_platform.Adaptive.simulation as si
from pele_platform.Utilities.BuildingBlocks import blocks as bb
import pele_platform.Utilities.Helpers.helpers as hp


class GPCRLauncher:

    def __init__(self, env):
        self.env = env

    def run(self):
        self.env.package = "gpcr_orth"

        result = bb.Pipeline([bb.GPCR, bb.Selection, bb.Rescoring], self.env).run()
        return result
