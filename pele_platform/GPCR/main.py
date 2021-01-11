from pele_platform.Utilities.Helpers.launcher_base import LauncherBase
from pele_platform.Utilities.BuildingBlocks.selection import Scatter6
from pele_platform.Utilities.BuildingBlocks.blocks import GPCR, Rescoring


class GPCRLauncher(LauncherBase):
    steps = [GPCR, Scatter6, Rescoring]
