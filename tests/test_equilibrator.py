import os

from pele_platform.main import run_platform_from_yaml
from .test_adaptive import test_path, check_file


def test_equilibrator_production():
    """
    Runs a full simulation with an additional short equilibration to estimate cluster values and conditions.
    """
    yaml_file = os.path.join(test_path, "induced_fit", "input_auto_clustering.yaml")
    job_params = run_platform_from_yaml(yaml_file)

    eq_output = os.path.join(job_params.pele_dir, "output", "preequilibration")
    assert os.path.isdir(eq_output)

    results = os.path.join(job_params.pele_dir, "results", "clusters")
    assert os.path.isdir(results)

    errors = check_file(job_params.pele_dir, "adaptive.conf", ['[1.48, 0.48]', '[2.0, 5, 7]'], [])
    assert not errors
