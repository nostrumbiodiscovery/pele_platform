from pele_platform.Utilities.Helpers.launcher_base import LauncherBase
from pele_platform.Utilities.BuildingBlocks.selection import LowestEnergy5
from pele_platform.Utilities.BuildingBlocks.simulation import InducedFitExhaustive, InducedFitFast, Rescoring


class InducedFitExhaustiveLauncher(LauncherBase):
    steps = [InducedFitExhaustive, LowestEnergy5, Rescoring]


class InducedFitFastLauncher(LauncherBase):
    steps = [InducedFitFast, LowestEnergy5, Rescoring]
