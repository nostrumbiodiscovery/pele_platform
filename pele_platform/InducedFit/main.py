from pele_platform.Utilities.Helpers.launcher_base import LauncherBase


class InducedFitExhaustiveLauncher(LauncherBase):
    steps = [{'type': 'InducedFitExhaustive'},
             {'type': 'LowestEnergy5'},
             {'type': 'Rescoring'}]


class InducedFitFastLauncher(LauncherBase):
    steps = [{'type': 'InducedFitFast'},
             {'type': 'LowestEnergy5'},
             {'type': 'Rescoring'}]
