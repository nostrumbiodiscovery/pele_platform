from pele_platform.Utilities.BuildingBlocks.pipeline import Pipeline


class LauncherBase:
    steps = []
    refinement_steps = []

    def __init__(self, env):
        self.env = env

    def run(self):
        steps = list(self.steps)  # Copy to prevent adding refinement steps more than once
        if not self.env.initial_args.skip_refinement:
            steps += self.refinement_steps

        result = Pipeline(steps, self.env).run()
        return result

