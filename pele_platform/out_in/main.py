from pele_platform.Utilities.Helpers.launcher_base import LauncherBase
from pele_platform.Utilities.BuildingBlocks.simulation import OutIn, Rescoring
from pele_platform.Utilities.BuildingBlocks.selection import ScatterN


class OutInLauncher(LauncherBase):
    steps = [OutIn, ScatterN, Rescoring]
