from pele_platform.Utilities.Helpers.launcher_base import LauncherBase
from pele_platform.Utilities.BuildingBlocks.selection import GMM
from pele_platform.Utilities.BuildingBlocks.blocks import Rescoring, InducedFitExhaustive


class PPILauncher(LauncherBase):
    steps = [InducedFitExhaustive]
    refinement_steps = [GMM, Rescoring]