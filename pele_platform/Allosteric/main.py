from pele_platform.Utilities.Helpers.launcher_base import LauncherBase
from pele_platform.Utilities.BuildingBlocks.selection import ScatterN
from pele_platform.Utilities.BuildingBlocks.simulation import GlobalExploration, InducedFitExhaustive, Rescoring


class AllostericLauncher(LauncherBase):
    steps = [GlobalExploration]
    refinement_steps = [ScatterN, InducedFitExhaustive, ScatterN, Rescoring]
