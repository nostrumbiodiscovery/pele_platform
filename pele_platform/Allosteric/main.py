from pele_platform.Utilities.Helpers.launcher_base import LauncherBase
from pele_platform.Utilities.BuildingBlocks.selection import Scatter6
from pele_platform.Utilities.BuildingBlocks.blocks import GlobalExploration, InducedFitExhaustive


class AllostericLauncher(LauncherBase):
    steps = [GlobalExploration]
    refinement_steps = [Scatter6, InducedFitExhaustive]