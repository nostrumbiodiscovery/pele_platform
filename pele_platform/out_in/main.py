from pele_platform.Utilities.Helpers.launcher_base import LauncherBase


class OutInLauncher(LauncherBase):
    steps = [{'type': 'OutIn'},
             {'type': 'ScatterN', 'options': {'distance': 6.0}},
             {{'type': 'Rescoring'}}]
