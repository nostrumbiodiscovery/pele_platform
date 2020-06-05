from . import test_adaptive as tk
import pele_platform.main as main
import pele_platform.Errors.custom_errors as ce
import pele_platform.constants.constants as cs
import os

test_path = os.path.join(cs.DIR, "Examples")
METAL_CONSTR_ARGS = os.path.join(test_path, "constraints/input_metals.yaml")
NO_METAL_CONSTR_ARGS = os.path.join(test_path, "constraints/input_no_metal_constraints.yaml")
FAIL_PERMISSIVE_METAL_CONSTR_ARGS = os.path.join(test_path, "constraints/fail_input_permissive_constraints.yaml")
PASS_PERMISSIVE_METAL_CONSTR_ARGS = os.path.join(test_path, "constraints/pass_input_permissive_constraints.yaml")
ALL_METAL_CONSTR_ARGS = os.path.join(test_path, "constraints/input_all_metal_constraints.yaml")
PERMISSIVE_EXCEPTION = os.path.join(test_path, "constraints/input_permissive_exception.yaml")

PASS_METAL_CONSTR = [
        '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.7238008975982666, "constrainThisAtom":  "A:40:_OG_", "toThisOtherAtom": "A:2007:MG__"}',
        '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 1.9963840246200562, "constrainThisAtom":  "Z:2001:_O5_", "toThisOtherAtom": "A:2007:MG__"},',
        '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.107039213180542, "constrainThisAtom":  "Z:2001:_O1_", "toThisOtherAtom": "A:2007:MG__"},',
        '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.232748031616211, "constrainThisAtom":  "A:17:_OG1", "toThisOtherAtom": "A:2007:MG__"},'
]

METAL_CONSTR = [
       '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.0319981575012207, "constrainThisAtom":  "A:239:_OD1", "toThisOtherAtom": "A:350:MG__"},',
       '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.0453453063964844, "constrainThisAtom":  "A:311:_OW_", "toThisOtherAtom": "A:350:MG__"},',
       '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.5390141010284424, "constrainThisAtom":  "A:401:CL__", "toThisOtherAtom": "A:350:MG__"},',
       '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.0555760860443115, "constrainThisAtom":  "A:312:_OW_", "toThisOtherAtom": "A:350:MG__"},',
       '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.0924105644226074, "constrainThisAtom":  "A:141:_OG_", "toThisOtherAtom": "A:350:MG__"},',
       '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.0959291458129883, "constrainThisAtom":  "A:139:_OG_", "toThisOtherAtom": "A:350:MG__"},'
]

ALL_METAL_CONSTR = [
        '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 1.920946717262268, "constrainThisAtom":  "A:268:_NE2", "toThisOtherAtom": "A:511:ZN__"},',
        '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.2836170196533203, "constrainThisAtom":  "A:609:_OW_", "toThisOtherAtom": "A:512:ZN__"},',
        '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.447920799255371, "constrainThisAtom":  "A:435:_NE2", "toThisOtherAtom": "A:512:ZN__"},',
        '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.749384641647339, "constrainThisAtom":  "A:766:_OW_", "toThisOtherAtom": "A:512:ZN__"},',
        '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.0511677265167236, "constrainThisAtom":  "A:294:_OE1", "toThisOtherAtom": "A:511:ZN__"},'
]


def test_metal_constraints(ext_args=METAL_CONSTR_ARGS):
    # checks metal constraints without any flags
    errors = []
    job, _ = main.run_platform(ext_args)
    errors = tk.check_file(job.pele_dir, "pele.conf", METAL_CONSTR, errors)
    assert not errors


def test_no_metal_constraints(ext_args=NO_METAL_CONSTR_ARGS):
    # checks no_metal_constraints flag
    errors = []
    job = main.run_platform(ext_args)
    errors = tk.check_file(job.pele_dir, "pele.conf", METAL_CONSTR, errors)
    assert errors 


def test_permissive_constraints(passed=PASS_PERMISSIVE_METAL_CONSTR_ARGS, failed=FAIL_PERMISSIVE_METAL_CONSTR_ARGS):
    
    # non-permissive -> supposed to fail due to lack of geometry
    try:
        job = main.run_platform(failed)
    except ce.NoGeometryAroundMetal:
        assert ce.NoGeometryAroundMetal
        return
    assert False

    # same system, but permissive -> should add constraints around the metal
    errors = []
    job = main.run_platform(passed)
    errors = tk.check_file(job.pele_dir, "pele.conf", PASS_METAL_CONSTR, errors)
    assert not errors


def test_all_metal_constraints(ext_args=ALL_METAL_CONSTR_ARGS, ext_args_permissive=PERMISSIVE_EXCEPTION):

    # checks constrain_all_metals -> should add whatever atoms in range
    errors = []
    job = main.run_platform(ext_args)
    errors = tk.check_file(job.pele_dir, "pele.conf", ALL_METAL_CONSTR, errors)
    assert not errors

    # same system, but permissive -> should fail due to lack of geometry
    try:
        job = main.run_platform(ext_args_permissive)
    except ce.NoGeometryAroundMetal:
        assert ce.NoGeometryAroundMetal
        return
    assert False

