import os

from pele_platform.Adaptive.equilibrator import Equilibrator
from pele_platform.main import run_platform_from_yaml
from .test_adaptive import test_path, check_file


def test_equilibrator_production():
    """
    Runs a full simulation with an additional short equilibration to estimate cluster values and conditions.
    """
    yaml_file = os.path.join(test_path, "induced_fit", "input_auto_clustering.yaml")
    job_params = run_platform_from_yaml(yaml_file)

    results = os.path.join(job_params.pele_dir, "results", "clusters")
    assert os.path.isdir(results)

    errors = check_file(job_params.pele_dir, "adaptive.conf", ['[2.0, 2.0, 2.0]', '[1.48, 1.48, 1.48]'], [])
    assert not errors


def test_values_extraction():
    """
    Checks extraction of cluster values and conditions given the path to the output folder.
    """
    output_path = os.path.join(
        test_path,
        "analysis",
        "data",
        "output",
    )

    cluster_values, cluster_conditions = Equilibrator.extract_cluster_values(
        output_path
    )

    assert cluster_conditions == [0.44, 0.244, 0.09]
    assert cluster_values == [2.5, 3.7, 4.0]
