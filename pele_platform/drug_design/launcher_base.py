from pele_platform.building_blocks.pipeline import Pipeline


class LauncherBase:
    # TODO: this would become Pipeline I guess? Confusing name
    steps = []
    refinement_steps = []

    def __init__(self, builder):
        self.builder = builder

    def run(self):
        steps = list(self.steps)  # Copy to prevent adding refinement steps more than once
        if not self.builder.initial_args.skip_refinement:
            steps += self.refinement_steps

        result = Pipeline(steps, self.builder).run()
        return result
