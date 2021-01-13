from pele_platform.Utilities.Helpers.launcher_base import LauncherBase


class PPILauncher(LauncherBase):
    steps = [{'type': 'InducedFitExhaustive'}]
    refinement_steps = [
        {'type': 'GMM'},
        {'type': 'Rescoring'}]
