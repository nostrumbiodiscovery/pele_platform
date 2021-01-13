from pele_platform.Utilities.Helpers.launcher_base import LauncherBase


class GPCRLauncher(LauncherBase):

    steps = [
        {'type': 'GPCR'},
        {'type': 'ScatterN', 'options': {"distance": 6.0}},
        {'type': 'Rescoring'}]
