import os
import pytest

from pele_platform import main
from pele_platform.constants import constants
from . import test_adaptive as ta

test_path = os.path.join(constants.DIR, "Examples")

expected_params1 = [
    '"type": "localNonBondingEnergy"',
    '"radius": 5',
    '"sideChainsToPerturb": { "links": {"ids": ["A:4"] } },',
    '"overlapFactor": 0.77,',
    '"numberOfStericTrials": 20,',
    '"numberOfTrials": 12,',
    '"gridResolution": 30,',
    '"atLeastOneSelectedTrial": true,'
]

expected_params2 = ['"refinementDistance": 15']


def test_covalent_docking_params():
    """
    Runs covalent docking in debug mode to ensure all metrics and user-defined variables are set correctly in pele.conf.
    """
    yaml_file = os.path.join(test_path, "covalent_docking", "input2.yaml")
    job, job2 = main.run_platform_from_yaml(yaml_file)

    errors = ta.check_file(job.pele_dir, "pele.conf", expected_params1, [])
    errors = ta.check_file(job2.pele_dir, "pele.conf", expected_params1+expected_params2, errors)
    assert not errors
