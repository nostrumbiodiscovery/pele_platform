import pytest

import pele_platform.Utilities.Parameters.pele_env as pv
import pele_platform.Errors.custom_errors as ce
from pele_platform.Utilities.BuildingBlocks.blocks import GPCR, Rescoring
from pele_platform.Utilities.BuildingBlocks.pipeline import Pipeline
from pele_platform.Utilities.BuildingBlocks.selection import Scatter6


@pytest.fixture
def pele_env():
    pele_env = pv.EnviroBuilder()
    return pele_env


@pytest.mark.parametrize("iterable", [([Scatter6, Rescoring]), ([GPCR, Rescoring]), ([])])
def test_pipeline_checker(pele_env, iterable):
    try:
        simulation_params = Pipeline(iterable, pele_env).run()
    except ce.PipelineError:
        assert True
        return
    assert False

