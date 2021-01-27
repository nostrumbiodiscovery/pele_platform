import os
import pytest
import pele_platform.constants.constants as cs
import pele_platform.main as main
from . import test_adaptive as ta

test_path = os.path.join(cs.DIR, "Examples")
BOX_RES_ARGS = os.path.join(test_path, "box/input.yaml")
GPCR_ARGS = os.path.join(test_path, "box/input_gpcr.yaml")

BOX_RES = [
    '"fixedCenter": [-3.6559998989105225, 64.46499633789062, 2.3980000019073486]',
]

GPCR_VALUES = [
    '"radius": 19.970223159033843,',
    '"fixedCenter": [-71.78435134887695, -13.431749963760375, -42.46209926605225]',
]


@pytest.mark.parametrize(
    ("yaml", "expected"), [(BOX_RES_ARGS, BOX_RES), (GPCR_ARGS, GPCR_VALUES)]
)
def test_box_center_from_residue(yaml, expected):
    """
    Runs simulation in debug mode to ensure the platform correctly identifies box center and radius in pele.conf.
    """
    errors = []
    (job,) = main.run_platform(yaml)
    errors = ta.check_file(job.pele_dir, "pele.conf", expected, errors)
    assert not errors
