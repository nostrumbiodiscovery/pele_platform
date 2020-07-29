import os
import pele_platform.constants.constants as cs
import pele_platform.main as main
from . import test_adaptive as ta

test_path = os.path.join(cs.DIR, "Examples")
BOX_RES_ARGS = os.path.join(test_path, "box/input.yaml")


BOX_RES = [
   '"fixedCenter": [-3.6559998989105225, 64.46499633789062, 2.3980000019073486]', 
]

def test_box_center_from_residue(ext_args=BOX_RES_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = ta.check_file(job.pele_dir, "pele.conf", BOX_RES, errors)
    assert not errors
