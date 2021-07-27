from pele_platform.building_blocks.pipeline import Pipeline


class LauncherBase:
    # TODO: this would become Pipeline I guess? Confusing name
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
