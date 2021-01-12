import inspect

import pele_platform.Adaptive.simulation as si
import pele_platform.features.adaptive as ft
from pele_platform.Utilities.Helpers.launcher_base import LauncherBase
import pele_platform.Utilities.BuildingBlocks.selection as selection
from pele_platform.Utilities.BuildingBlocks.selection import *
from pele_platform.Utilities.BuildingBlocks.blocks import *
import pele_platform.Utilities.BuildingBlocks.blocks as blocks
import pele_platform.Errors.custom_errors as ce


class AdaptiveLauncher(LauncherBase):

    def run(self):
        self.set_simulation_type()
        self.env.build_adaptive_variables(self.env.initial_args)
        self.env.create_files_and_folders()
        self.env = si.run_adaptive(self.env)
        return self.env

    def set_simulation_type(self):
        # NEEDS IMPROVEMENT. Ensuring it doesn't crash in features.adaptive
        # with something like ('EnviroBuilder' object has no attribute 'full'), do you have a better idea?
        for arg, value in self.env.initial_args.__dict__.items():
            if arg in ft.all_simulations:
                setattr(self.env, arg, value)


class WorkflowLauncher(LauncherBase):

    @property
    def steps(self):
        available = {**dict((name, func) for name, func in inspect.getmembers(selection)), **dict((name, func) for name, func in inspect.getmembers(blocks))}
        iterable = self.env.initial_args.workflow
        for i in iterable:
            if not (i in available.keys() or inspect.isclass(i)):
                raise ce.PipelineError(
                    "Block {} cannot be found. Please check spelling and refer to the PELE Platform documentation "
                    "for an up-to-date list of available BuildingBlocks".format(i))

        return [eval(elem) for elem in iterable]
