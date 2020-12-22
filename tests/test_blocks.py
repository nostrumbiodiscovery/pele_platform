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
    '"type" : "epsilon"',
    '"iterations" : 50',
    '"peleSteps" : 8',
    '"conditions": [0.7, 0.4, 0.0]'
]

OutIn_lines = [
    '"type" : "inverselyProportional"',
    '"iterations" : 100',
    '"peleSteps" : 8',
    '"conditions": [1, 0.6, 0.0]'
]

GlobalExploration_lines = [
    '"type" : "inverselyProportional"',
    '"iterations" : 100',
    '"peleSteps" : 4',
    '"conditions": [1, 0.6, 0.0]'
]

InducedFitExhaustive_lines = [
    '"type" : "independent"',
    '"iterations" : 1',
    '"peleSteps" : 1000',
    '"conditions": [1, 0.6, 0.4, 0.0]'
]

InducedFitFast_lines = [
    '"type" : "inverselyProportional"',
    '"iterations" : 30',
    '"peleSteps" : 12',
    '"conditions": [1, 0.6, 0.4, 0.0]'
]

Rescoring_lines = [
    '"type" : "independent"',
    '"iterations" : 20',
    '"peleSteps" : 12',
    '"conditions": [1, 0.6, 0.4, 0.0]'
]


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
                         [("gpcr/input_defaults.yaml", "gpcr_orth", bb.GPCR, GPCR_lines),
                          ("global/input_defaults.yaml", "global", bb.GlobalExploration, GlobalExploration_lines),
                          ("out_in/input_default.yaml", "out_in", bb.OutIn, OutIn_lines),
                          ("induced_fit/input_exhaustive_defaults.yaml", "induced_fit_exhaustive", bb.InducedFitExhaustive,
                           InducedFitExhaustive_lines),
                          ("induced_fit/input_fast_defaults.yaml", "induced_fit_fast", bb.InducedFitFast,
                           InducedFitFast_lines),
                          ("rescoring/input_defaults.yaml", "rescoring", bb.Rescoring, Rescoring_lines)])
def test_simulation_blocks(yaml, package, block, expected):
    
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
    params = simulation_block.run()
    
    directory = params.pele_dir
    errors = []
    errors = tk.check_file(directory, "adaptive.conf", expected, errors)
    assert not errors

