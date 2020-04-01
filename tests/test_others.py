import os
import glob
import pele_platform.constants.constants as cs
import pele_platform.constants.pele_params as pp
import pele_platform.main as main
import test_adaptive as tk
import pytest

test_path = os.path.join(cs.DIR, "Examples")
EXTERNAL_CONSTR_ARGS = os.path.join(test_path, "constraints/input_external_constraints.yaml")
PPP_CONSTR_ARGS = os.path.join(test_path, "constraints/input_ppp.yaml")


EXT_CONSTR = [
    '{ "type": "constrainAtomToPosition", "springConstant": 10, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_H__" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_H__" },',
    '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.34, "constrainThisAtom":  "A:1:_H__", "toThisOtherAtom": "L:1:_C21"},',
    '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.34, "constrainThisAtom":  "A:1:_H__", "toThisOtherAtom": "L:1:_C21"}'
]

PPP_CONSTR = [
    '"constrainThisAtom": "B:207:_CA_" }',
    '"constrainThisAtom": "A:1:_CA_" }',
    '"constrainThisAtom": "B:247:_CA_" }'
]

def test_external_constraints(ext_args=EXTERNAL_CONSTR_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = tk.check_file(job.pele_dir, "pele.conf", EXT_CONSTR, errors)
    assert not errors

def test_ppp_constraints(ext_args=PPP_CONSTR_ARGS):
    errors = []
    job = main.run_platform(ext_args)
    errors = tk.check_file(job.pele_dir, "pele.conf", PPP_CONSTR, errors)
    assert not errors

def test_checker():
    yaml = os.path.join(test_path, "checker/input.yaml")
    try:
        job = main.run_platform(yaml)
    except KeyError:
        assert KeyError
        return
    assert False
