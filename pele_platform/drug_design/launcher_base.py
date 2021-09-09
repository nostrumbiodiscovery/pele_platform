from pele_platform.building_blocks.pipeline import Pipeline
from pele_platform.context import context


class LauncherBase:
    # TODO: this would become Pipeline I guess? Confusing name
    steps = []
    refinement_steps = []

    def run(self):
        steps = list(self.steps)  # Copy to prevent adding refinement steps more than once
        if not context.yaml_parser.skip_refinement:
            steps += self.refinement_steps

        result = Pipeline(steps).run()
        return result
