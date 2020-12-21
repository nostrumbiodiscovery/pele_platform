import os
import pytest

import pele_platform.Utilities.Parameters.pele_env as pv
import pele_platform.Errors.custom_errors as ce
from pele_platform.constants import constants as cs
import pele_platform.Utilities.BuildingBlocks.blocks as bb
from pele_platform.Utilities.BuildingBlocks.pipeline import Pipeline
from pele_platform.Utilities.BuildingBlocks.selection import Scatter6
import pele_platform.Checker.valid_flags as vf
import pele_platform.Utilities.Helpers.yaml_parser as yp
from . import test_adaptive as tk

test_path = os.path.join(cs.DIR, "Examples")

GPCR_lines = [
    '"spawning_type": "epsilon"',
    '"iterations": 50',
    '"pele_steps": 8',
    '"cluster_conditions": "[0.7, 0.4, 0.0]"'
]

OutIn_lines = [
    '"spawning_type": "inverselyProportional"',
    '"iterations": 100',
    '"pele_steps": 8',
    '"cluster_conditions": "[1, 0.6, 0.0]"'
]

GlobalExploration_lines = [
    '"spawning_type": "inverselyProportional"',
    '"iterations": 100',
    '"pele_steps": 4',
    '"cluster_conditions": "[1, 0.6, 0.0]"'
]

InducedFitExhaustive_lines = [
    '"spawning_type": "independent"',
    '"iterations": 1',
    '"pele_steps": 1000',
    '"cluster_conditions": "[1, 0.6, 0.4, 0.0]"'
]

InducedFitFast_lines = [
    '"spawning_type": "inverselyProportional"',
    '"iterations": 30',
    '"pele_steps": 12',
    '"cluster_conditions": "[1, 0.6, 0.4, 0.0]"'
]

Rescoring_lines = [
    '"spawning_type": "independent"',
    '"iterations": 20',
    '"pele_steps": 12',
    '"cluster_conditions": "[1, 0.6, 0.4, 0.0]"'
]


@pytest.fixture
def pele_env(yaml, package, block):
    # get YamlParser ready
    yaml_file = os.path.join(test_path, yaml)
    yaml = yp.YamlParser(yaml_file, vf.VALID_FLAGS_PLATFORM)
    yaml.read()

    # create pele environment
    pele_env = pv.EnviroBuilder()
    pele_env.initial_args = yaml
    pele_env.initial_args.package = pele_env.package = package

    # run Building Block
    simulation_block = block(pele_env, "test_folder")
    simulation_params = simulation_block.run()

    return simulation_params


@pytest.mark.parametrize("iterable", [([Scatter6, bb.Rescoring]), ([bb.GPCR, bb.Rescoring]), ([])])
def test_pipeline_checker(iterable):
    try:
        env = pv.EnviroBuilder()
        simulation_params = Pipeline(iterable, env).run()
    except ce.PipelineError:
        assert True
        return
    assert False


@pytest.mark.parametrize(("yaml", "package", "block", "expected"),
                         [("gpcr/input.yaml", "gpcr_orth", bb.GPCR, GPCR_lines),
                          ("global/input.yaml", "global", bb.GlobalExploration, GlobalExploration_lines),
                          ("out_in/input.yaml", "out_in", bb.OutIn, OutIn_lines),
                          ("induced_fit/input_exhaustive.yaml", "induced_fit_exhaustive", bb.InducedFitExhaustive,
                           InducedFitExhaustive_lines), (
                                  "induced_fit/input_fast.yaml", "induced_fit_fast", bb.InducedFitFast,
                                  InducedFitFast_lines),
                          ("rescoring/input.yaml", "rescoring", bb.Rescoring, Rescoring_lines)])
def test_simulation_blocks(pele_env, yaml, package, block, expected):
    params = pele_env(yaml, package, block)
    directory = params.pele_dir
    errors = []
    errors = tk.check_file(directory, "pele.conf", expected, errors)
    assert not errors
