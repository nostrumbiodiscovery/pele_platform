import os
import glob
import pele_platform.constants.constants as cs
import pele_platform.constants.pele_params as pp
import pele_platform.main as main
import test_kinase as tk

test_path = os.path.join(cs.DIR, "Examples")
EXTERNAL_CONSTR_ARGS = [os.path.join(test_path, "induced_fit/input_external_constraints.yaml")]


EXT_CONSTR = [
    '{ "type": "constrainAtomToPosition", "springConstant": 10, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_H__" },',
    '{ "type": "constrainAtomToPosition", "springConstant": 5, "equilibriumDistance": 0.0, "constrainThisAtom": "A:1:_H__" },',
    '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.34, "constrainThisAtom":  "A:1:_H__", "toThisOtherAtom": "L:1:_C21"},',
    '{"type": "constrainAtomsDistance", "springConstant": 50, "equilibriumDistance": 2.34, "constrainThisAtom":  "A:1:_H__", "toThisOtherAtom": "L:1:_C21"}'
]
def test_external_constraints(ext_args=EXTERNAL_CONSTR_ARGS):
    errors = []
    arguments = main.parseargs_yaml(ext_args)
    arguments = main.YamlParser(arguments.input_file)
    job = main.Launcher(arguments).launch()
    errors = tk.check_file(job.pele_dir, "pele.conf", EXT_CONSTR, errors)
    assert not errors

