from pele_platform.Utilities.Helpers.launcher_base import LauncherBase


class AllostericLauncher(LauncherBase):
    steps = [{'type': 'GlobalExploration'}]

    refinement_steps = [
        {'type': 'ScatterN', 'options': {"distance": 6.0}},
        {'type': 'InducedFitExhaustive'},
        {'type': 'ScatterN', 'options': {"distance": 3.0}},
        {'type': 'Rescoring'}]
