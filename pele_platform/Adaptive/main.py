import pele_platform.Adaptive.simulation as si
import pele_platform.features.adaptive as ft


class AdaptiveLauncher:

    def __init__(self, env):
        self.env = env
        self.env.package = "adaptive"

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
