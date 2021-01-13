from pele_platform.Utilities.Helpers.launcher_base import LauncherBase
from pele_platform.Utilities.BuildingBlocks.selection import ScatterN
from pele_platform.Utilities.BuildingBlocks.simulation import GPCR, Rescoring


class GPCRLauncher(LauncherBase):
    steps = [GPCR, ScatterN, Rescoring]
