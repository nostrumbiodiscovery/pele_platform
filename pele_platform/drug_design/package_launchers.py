import inspect

from pele_platform.Adaptive import simulation as si
from pele_platform.drug_design.launcher_base import LauncherBase
from pele_platform.Errors import custom_errors as ce
from pele_platform.building_blocks import simulation as blocks
from pele_platform.building_blocks import selection as selection
from pele_platform.features import adaptive as ft
from pele_platform.context import context


class AdaptiveLauncher(LauncherBase):

    def run(self):
        self.set_simulation_type()
        context.parameters_builder.build_adaptive_variables()
        context.parameters.create_files_and_folders()
        si.run_adaptive()
        return [context.parameters]

    def set_simulation_type(self):
        # NEEDS IMPROVEMENT. Ensuring it doesn't crash in features.adaptive
        # with something like ('EnviroBuilder' object has no attribute 'full'), do you have a better idea?
        for arg in dir(context.yaml_parser):
            if arg in ft.all_simulations:
                setattr(context.parameters_builder, arg, getattr(context.yaml_parser, arg))


class WorkflowLauncher(LauncherBase):

    @property
    def steps(self):

        available = {**dict((name, func) for name, func in inspect.getmembers(selection)),
                     **dict((name, func) for name, func in inspect.getmembers(blocks))}

        iterable = context.yaml_parser.workflow
        simulation_blocks = [i.get('type', None) for i in iterable]

        for i in simulation_blocks:
            if not (i in available.keys() or inspect.isclass(i)):
                raise ce.PipelineError(
                    "Block {} cannot be found. Please check spelling and refer to the PELE Platform documentation "
                    "for an up-to-date list of available building_blocks".format(i))
        return iterable


class SiteFinderLauncher(LauncherBase):
    steps = [{'type': 'GlobalExploration'}]
    refinement_steps = [
        {'type': 'Clusters'},
        {'type': 'LocalExplorationExhaustive'},
        {'type': 'Clusters'},
        {'type': 'Rescoring'}]


class GPCRLauncher(LauncherBase):
    steps = [
        {'type': 'GPCR'},
        {'type': 'Clusters'},
        {'type': 'Rescoring'}]


class InducedFitFastLauncher(LauncherBase):
    steps = [{'type': 'LocalExplorationFast'},
             {'type': 'LowestEnergy'},
             {'type': 'Rescoring'}]


class InducedFitExhaustiveLauncher(LauncherBase):
    steps = [{'type': 'LocalExplorationExhaustive'},
             {'type': 'LowestEnergy'},
             {'type': 'Rescoring'}]


class OutInLauncher(LauncherBase):
    steps = [{'type': 'OutIn'},
             {'type': 'Clusters'},
             {'type': 'Rescoring'}]


class PPILauncher(LauncherBase):
    steps = [{'type': 'LocalExplorationExhaustive'}]
    refinement_steps = [
        {'type': 'GMM'},
        {'type': 'Rescoring'}]


class CovalentDocking(LauncherBase):
    steps = [{'type': 'CovalentDockingExploration'}]
    refinement_steps = [
        {'type': 'LowestLocalNonbondingEnergy'},
        {'type': 'CovalentDockingRefinement'}]
