import os
import pytest

from pele_platform.main import run_platform_from_yaml
from .test_adaptive import test_path


@pytest.mark.parametrize(
    ("yaml_file", "conditions", "values"),
    [
        (
            os.path.join("induced_fit", "input_auto_clustering.yaml"),
            [1.48, 0.48],
            [2.0, 5, 7],
        ),
        (
            os.path.join("out_in", "input_auto_clustering.yaml"),
            [0.44, 0.0],
            [2, 5, 7],
        ),
    ],
)
def test_equilibrator_production(yaml_file, conditions, values):
    """
    Runs a full simulation with an additional short equilibration to estimate cluster values and conditions.
    """
    yaml_file = os.path.join(test_path, yaml_file)
    job_params = run_platform_from_yaml(yaml_file)

    eq_output = os.path.join(job_params.pele_dir, "output", "preequilibration")
    assert os.path.isdir(eq_output)

    results = os.path.join(job_params.pele_dir, "results", "clusters")
    assert os.path.isdir(results)

    assert job_params.cluster_conditions == conditions
    assert job_params.cluster_values == values
