import pele_platform.Utilities.Parameters.pele_env as pv
import pele_platform.Errors.custom_errors as ce
from pele_platform.Utilities.BuildingBlocks.blocks import GPCR, Rescoring
from pele_platform.Utilities.BuildingBlocks.pipeline import Pipeline
from pele_platform.Utilities.BuildingBlocks.selection import Scatter6


pele_env = pv.EnviroBuilder()


def test_pipeline_checker():

    try:
        pipeline_wrong_start = Pipeline([Scatter6, Rescoring], pele_env).run()
    except ce.PipelineError as e:
        assert str(e) == "Pipeline should begin and end with a Simulation block, e.g. GlobalExploration."
        return
    assert False, "Test did not catch Pipeline starting with Selection block."


def test_pipeline_checker_2():

    try:
        pipeline_no_selection = Pipeline([GPCR, Rescoring], pele_env).run()
    except ce.PipelineError as e:
        assert True
        return
    assert False, "Test did not catch Pipeline without Selection blocks."


def test_empty_pipeline():

    try:
        pipeline_empty = Pipeline([], pele_env).run()
    except ce.PipelineError as e:
        assert str(e).strip("'") == "Pipeline doesn't contain any BuildingBlocks."
        return
    assert False, "Test did not catch empty Pipeline error."


