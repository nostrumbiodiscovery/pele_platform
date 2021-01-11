from pele_platform.Utilities.Helpers.launcher_base import LauncherBase
from pele_platform.Utilities.BuildingBlocks.blocks import OutIn, Rescoring
from pele_platform.Utilities.BuildingBlocks.selection import Scatter6


class OutInLauncher(LauncherBase):
    steps = [OutIn, Scatter6, Rescoring]
